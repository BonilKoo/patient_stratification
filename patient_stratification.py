import argparse
import os
import random
import copy
import numpy as np
import pandas as pd
import parmap
from scipy.stats import kendalltau
from sklearn.linear_model import LogisticRegression
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceError
import warnings
warnings.filterwarnings(action='ignore', category=RuntimeWarning)

def parse_args():
    parser = argparse.ArgumentParser()
    
    # - -- type required default help
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--expression', required=True, help='gene x sample')
    parser.add_argument('--subtype', required=True, help='sample subtype (no header)')
    parser.add_argument('--survival', required=True, help='sample event duration (no header)')
    parser.add_argument('--subtype_order', required=True, type=str, help='sparated by comma')
    parser.add_argument('--survival_coef', type=float, default=0.5)
    parser.add_argument('--n_threads', type=int, default=1)
    parser.add_argument('--n_subpop', type=int, default=10)
    parser.add_argument('--n_chrom', type=int, default=150)
    parser.add_argument('--n_init_gene', type=int, default=10)
    parser.add_argument('--elite_proportion', type=float, default=0.1)
    parser.add_argument('--num_participant', type=int, default=16)
    parser.add_argument('--crossover_prob', type=float, default=0.3)
    parser.add_argument('--mutation_prob', type=float, default=0.0005)
    parser.add_argument('--patience', type=int, default=10)

    return parser.parse_args()

def set_seed(seed):
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)

def rmfile(filename):
    if os.path.exists(filename):
        os.system(f'rm -rf {filename}')

def mkdir(dir_name):
    os.system(f'mkdir -p {dir_name}')

def load_expr_data(expr_file):
    expr = pd.read_table(expr_file, index_col=0)
    expr_geneset = list(expr.index)
    return expr, expr_geneset

def load_subtype_data(subtype_file):
    return pd.read_table(subtype_file, names=['sample', 'subtype'])

def load_survival_data(survival_file):
    return pd.read_table(survival_file, names=['sample', 'event', 'duration'])

class chromosome:
    def __init__(self, n_init_gene):
        selected_genes = np.random.choice(globals()['geneset'], size=n_init_gene, replace=False)
        ascending = np.random.choice([True, False], size=n_init_gene, replace=True)
        self.chrom_info = {selected_genes[idx]:ascending[idx] for idx in range(n_init_gene)}

    def sample_rank_in_chrom(self):
        sample_rank = pd.DataFrame([globals()['expr'].loc[gene].rank(ascending=self.chrom_info[gene]) for gene in self.chrom_info.keys()])
        sample_rank = pd.DataFrame(sample_rank.mean().sort_values())
        sample_rank.columns = ['risk_score']
        sample_rank = sample_rank.merge(globals()['subtype_info'], left_index=True, right_on='sample')
        return sample_rank

    def kendall_tau_b(self, sample_rank):
        target_subtype_list = []
        num = 1
        for subtype in globals()['subtype_order']:
            target_subtype_list.extend([num] * sample_rank['subtype'].value_counts()[subtype])
            sample_rank.loc[sample_rank['subtype'] == subtype, 'subtype_num'] = num
            num += 1
        return kendalltau(target_subtype_list, sample_rank['subtype_num']).correlation

    def set_boundary(self, sample_rank):
        classifier = LogisticRegression(penalty='none', class_weight='balanced').fit(sample_rank['risk_score'].values.reshape(-1,1), sample_rank['subtype'])
        predicted = classifier.predict(sample_rank['risk_score'].values.reshape(-1,1))
    
        try:
            boundary_list = [(np.where(predicted == subtype)[0][0], np.where(predicted == subtype)[0][-1]) for subtype in globals()['subtype_order']]
        except IndexError:
            boundary_list = None

        return boundary_list

    def devide_series(self, subtype, boundary, sample_rank):
        intermediate = None
        if boundary[0] == 0:
            low = sample_rank.iloc[:boundary[1]+1]
            high = sample_rank.iloc[boundary[1]+1:]
        elif boundary[1] == len(sample_rank)-1:
            low = sample_rank.iloc[:boundary[0]]
            high = sample_rank.iloc[boundary[0]:]
        else:
            low = sample_rank.iloc[:boundary[0]]
            intermediate = sample_rank.iloc[boundary[0]:boundary[1]+1]
            high = sample_rank.iloc[boundary[1]+1:]

        low = low[low['subtype'] == subtype]['sample'].to_list()
        if intermediate is not None:
            intermediate = intermediate[intermediate['subtype'] == subtype]['sample'].to_list()
        high = high[high['subtype'] == subtype]['sample'].to_list()

        return low, intermediate, high

    def cal_survival_score(self, survival_low, survival_high):
        logrank_p = logrank_test(survival_low['duration'], survival_high['duration'], survival_low['event'], survival_high['event']).p_value
        if logrank_p < 0.01:
            logrank_p = 0.01 # upperbound

        survival_low = survival_low.copy()
        survival_high = survival_high.copy()
        survival_low.loc[:, 'group'] = ['low'] * len(survival_low)
        survival_high.loc[:, 'group'] = ['high'] * len(survival_high)
        survival_group = survival_low.append(survival_high)

        cph = CoxPHFitter()
        try:
            cph.fit(survival_group, duration_col='duration', event_col='event', formula='group')
        except ConvergenceError:
            return -2
        cph_result = cph.summary

        if cph_result.index[0] == 'group[T.high]':
            if cph_result['coef'][0] >= 0:
                score = np.log10(logrank_p) * -1
            else:
                score = np.log10(logrank_p) * 1
        else:
            if cph_result['coef'][0] >= 0:
                score = np.log10(logrank_p) * 1
            else:
                score = np.log10(logrank_p) * -1

        return score

    def survival_score(self, low, intermediate, high):
        survival_low = globals()['survival_info'][globals()['survival_info']['sample'].isin(low)]
        if intermediate is not None:
            survival_intermediate = globals()['survival_info'][globals()['survival_info']['sample'].isin(intermediate)]
        survival_high = globals()['survival_info'][globals()['survival_info']['sample'].isin(high)]

        if intermediate is None:
            if (len(low) == 0) or (len(high) == 0):
                return [0]
            else:
                return [self.cal_survival_score(survival_low, survival_high)]

        elif ((len(low) == 0) and (len(intermediate) == 0)) or ((len(intermediate) == 0) and (len(high) == 0)) or ((len(low) == 0) and (len(high) == 0)):
            return [0, 0]

        elif len(low) == 0:
            return [self.cal_survival_score(survival_intermediate, survival_high)]
        elif len(intermediate) == 0:
            return [self.cal_survival_score(survival_low, survival_high)]
        elif len(high) == 0:
            return [self.cal_survival_score(survival_low, survival_intermediate)]

        else:
            return [self.cal_survival_score(survival_low, survival_intermediate),
                   self.cal_survival_score(survival_intermediate, survival_high)]

    def fitness_function(self, survival_flag=0):
        set_seed(globals()[f'seed']) #####
        if len(self.chrom_info) == 0: # check # selected nodes == 0
            return -1

        sample_rank = self.sample_rank_in_chrom()

        kendalltau_correlation = self.kendall_tau_b(sample_rank)

        if survival_flag == 0:
            return kendalltau_correlation

        boundary_list = self.set_boundary(sample_rank)
        if boundary_list is None:
            return -2

        survival_scores = []
        for idx, subtype in enumerate(globals()['subtype_order']):
            low, intermediate, high = self.devide_series(subtype, boundary_list[idx], sample_rank)

            score = self.survival_score(low, intermediate, high)
            survival_scores.extend(score)

        survival_score = sum(survival_scores) / len(survival_scores)

        return kendalltau_correlation + globals()['survival_coef'] * survival_score

    def mutation(self, mutation_prob):
        mutation_flag = np.random.choice([True, False], size=len(globals()['geneset']), p=[mutation_prob, 1-mutation_prob])
        mutated_genes = np.array(globals()['geneset'])[mutation_flag]
        for gene in mutated_genes:
            if gene in self.chrom_info.keys():
                del self.chrom_info[gene]
            else:
                self.chrom_info[gene] = np.random.choice([True, False], size=1)[0]

class subpopulation:
    def __init__(self, n_chrom, n_init_gene):
        self.chrom_list = [chromosome(n_init_gene) for _ in range(n_chrom)]

    def fitness_function(self, survival_flag):
        chrom_fitness_list = parmap.map(chromosome.fitness_function, self.chrom_list, survival_flag, pm_processes=globals()['n_threads'])
        self.chrom_fitness_dict = dict(zip(self.chrom_list, chrom_fitness_list))

    def elitism(self, ratio):
        chrom_fitness_df = pd.DataFrame(self.chrom_fitness_dict, index=['score']).T.sort_values(by='score', ascending=False)
        self.top = chrom_fitness_df.index[0]
        self.elites = chrom_fitness_df.index[:round(len(chrom_fitness_df) * ratio)].to_list()
        self.non_elites = chrom_fitness_df.index[round(len(chrom_fitness_df) * ratio):].to_list()
        self.bottom = chrom_fitness_df.index[-1]

    def tournament_selection(self, num_participant):
        set_seed(globals()[f'seed']) #####
        participants = np.random.choice(self.non_elites, size=num_participant, replace=False)
        participants_fitness_dict = {participant:self.chrom_fitness_dict[participant] for participant in participants}
        winner = list(participants_fitness_dict.keys())[list(participants_fitness_dict.values()).index(max(participants_fitness_dict.values()))]
        return winner

    def selection_in_non_elites(self, num_participant):
        return parmap.map(self.tournament_selection, len(self.non_elites) * [num_participant], pm_processes=globals()['n_threads'])

    def crossover(self, idx, selected_chroms, crossover_prob):
        set_seed(globals()[f'seed']) #####
        chrom1 = copy.deepcopy(selected_chroms[2*idx])
        chrom2 = copy.deepcopy(selected_chroms[2*idx+1])

        chrom1_genes = set(chrom1.chrom_info.keys())
        chrom2_genes = set(chrom2.chrom_info.keys())

        chrom1_only = list(chrom1_genes - chrom2_genes)
        chrom2_only = list(chrom2_genes - chrom1_genes)

        for gene in chrom1_only:
            if np.random.choice([True, False], size=1, p=[crossover_prob, 1-crossover_prob])[0]:
                ascending = chrom1.chrom_info[gene]
                del chrom1.chrom_info[gene]
                chrom2.chrom_info[gene] = ascending
        for gene in chrom2_only:
            if np.random.choice([True, False], size=1, p=[crossover_prob, 1-crossover_prob])[0]:
                ascending = chrom2.chrom_info[gene]
                del chrom2.chrom_info[gene]
                chrom1.chrom_info[gene] = ascending

        # chrom1 chrom2 intersection

        return [chrom1, chrom2]

    def crossover_in_non_elites(self, selected_chroms, crossover_prob):
        return parmap.map(self.crossover, range(int(len(selected_chroms)/2)), selected_chroms, crossover_prob, pm_processes=globals()['n_threads'])

    def selection_crossover_mutation(self, num_participant, crossover_prob, mutation_prob):
        next_generation = self.elites

        # deterministic tournament selection
        selected_chroms = self.selection_in_non_elites(num_participant)

        # uniform crossover
        if len(selected_chroms) % 2 == 0:
            crossovered_chroms = np.array(self.crossover_in_non_elites(selected_chroms, crossover_prob)).reshape(-1).tolist()
        else:
            crossovered_chroms = np.array(self.crossover_in_non_elites(selected_chroms[:-1], crossover_prob)).reshape(-1).tolist()
            crossovered_chroms.append(selected_chroms[-1])

        # bit flip mutation
        for chrom in crossovered_chroms:
            chrom.mutation(mutation_prob)

        next_generation += crossovered_chroms
        self.next_generation = next_generation

    def update_chrom_list(self):
        self.chrom_list = self.next_generation

class population:
    def __init__(self, n_subpop, n_chrom, n_init_gene):
        self.subpop_list = [subpopulation(n_chrom, n_init_gene) for _ in range(n_subpop)]

    def fitness_function(self, survival_flag):
        for subpop in self.subpop_list:
            subpop.fitness_function(survival_flag)

    def elitism(self, ratio, flag=0):
        self.ratio = ratio
        for subpop in self.subpop_list:
            subpop.elitism(ratio)

        subpop_top_fitness_dict  = {subpop.top:subpop.chrom_fitness_dict[subpop.top] for subpop in self.subpop_list}
        self.top_fitness = max(subpop_top_fitness_dict.values())
        self.top_chrom = list(subpop_top_fitness_dict.keys())[list(subpop_top_fitness_dict.values()).index(self.top_fitness)]

        if flag == 0:
            for subpop in self.subpop_list:
                subpop.bottom.chrom_info = self.top_chrom.chrom_info
                subpop.chrom_fitness_dict[subpop.bottom] = self.top_fitness

    def selection_crossover_mutation(self, num_participant, crossover_prob, mutation_prob):
        for subpop in self.subpop_list:
            subpop.selection_crossover_mutation(num_participant, crossover_prob, mutation_prob)

    def next_generation(self):
        for subpop in self.subpop_list:
            subpop.update_chrom_list()

    def check_improvement(self, patience):
        before_top_fitness = self.top_fitness

        self.elitism(self.ratio, flag=1)
        after_top_fitness = self.top_fitness

        if after_top_fitness > before_top_fitness:
            return 0
        else:
            return patience+1

    def save_best_chrom_info(self, filename):
        best_chrom = pd.DataFrame(self.top_chrom.chrom_info, index=['ascending']).T
        best_chrom.index.name = 'gene'
        best_chrom.to_csv(filename, sep='\t')

def save_sample_info(chrom_file, sample_file):
    chrom_info = pd.read_table(chrom_file)
    chrom_info = {chrom_info.loc[idx]['gene']:chrom_info.loc[idx]['ascending'] for idx in chrom_info.index}
    
    res = chromosome(0)
    res.chrom_info = chrom_info
    
    rank_table = res.sample_rank_in_chrom(expr, subtype_info).reset_index(drop=True)
    boundary_list = res.set_boundary(rank_table, flag=0)
    
    subtype_order = globals()[f'subtype_order']
    subtype_idx_dict = {subtype_order[idx]:idx for idx in range(len(subtype_order))}
    
    output_file = open(sample_file, 'w')
    output_file.write('sample\tsubtype\trisk\n')
    for idx in rank_table.index:
        row = rank_table.loc[idx]
        sample = row['sample']        
        subtype = row['subtype']
        subtype_idx = subtype_idx_dict[subtype]
        if (subtype_idx == 0) or (subtype_idx == len(subtype_order)-1):
            if boundary_list[0][0] <= idx <= boundary_list[subtype_idx-1][1]:
                risk = 'Low'
            else:
                risk = 'High'
        else:
            if boundary_list[0][0] <= idx <= boundary_list[subtype_indx-1][1]:
                risk = 'Low'
            elif boundary_list[subtype_idx][0] <= idx <= boundary_list[subtype_idx][1]:
                risk = 'Intermediate'
            else:
                risk = 'High'
        output_file.write(f'{sample}\t{subtype}\t{risk}')
    output_file.close()

def main():
    args = parse_args()

    set_seed(args.seed)
    rmfile(args.outdir)
    mkdir(args.outdir)
#     mkdir(args.outdir+'/best_geneset')

    globals()['expr'], globals()['geneset'] = load_expr_data(args.expression)
    globals()['subtype_info'] = load_subtype_data(args.subtype)
    globals()['survival_info'] = load_survival_data(args.survival)

    globals()['subtype_order'] = args.subtype_order.split(',')
    globals()['n_threads'] = args.n_threads
    globals()['survival_coef'] = args.survival_coef
    survival_flag = 0
    patience = 0
    iteration = 0
    
    print('Initialization...')
    pop = population(n_subpop=args.n_subpop, n_chrom=args.n_chrom, n_init_gene=args.n_init_gene)

    while True:
        iteration += 1
        print(f'Iteration {iteration}...')
        pop.fitness_function(survival_flag)
        pop.elitism(args.elite_proportion, flag=0)

#         pop.save_best_chrom_info(args.outdir+f'/best_geneset/result{iteration}.txt')

        pop.selection_crossover_mutation(args.num_participant, args.crossover_prob, args.mutation_prob)
        pop.next_generation()
        pop.fitness_function(survival_flag)
        patience = pop.check_improvement(patience)
        if patience == args.patience:
            break

        if (survival_flag == 0) and (((np.array([list(subpop.chrom_fitness_dict.values()) for subpop in pop.subpop_list]).reshape(-1) > 0.75).sum() / (args.n_subpop * args.n_chrom)) > 0.95) and (survival_coef != 0):
            survival_flag = 1
            patience = 0
            print('Start including survival term...')

    print('END...')
    best_geneset_file = args.outdir+'/geneset.txt'
    sample_result_file = args.outdir+'/sample_risk.txt'
    pop.save_best_chrom_info(best_geneset_file)
    save_sample_info(best_geneset_file, sample_result_file)

if __name__ == '__main__':
    args = parse_args()
    seed = args.seed
    main()