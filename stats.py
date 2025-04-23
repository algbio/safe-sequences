from collections import defaultdict
import argparse
import re


def parse_input_file(filename):
    data = defaultdict(list)
    current_graph = None

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('#Graph'):
                current_graph = int(line.split()[1])

            elif current_graph is not None:
                if re.match(r'\d+,\s*\d+,\s*\d+', line):
                    n, m, w = map(int, line.split(','))

                    data[current_graph] = {
                        'n': n, 'm': m, 'w': w,
                        'solved default': None, 'total time default': None,
                        'solved sequences heur': None, 'total time sequences heur': None,
                        'preprocess sequences heur': None, 'ilp time seqs heur': None,
                        'fixed vars seqs': None,
                    }
                
                elif 'solved' in line or 'time' in line or 'preprocess' in line or 'fixed' in line:

                    key, value = line.split(':')
                    key = key.strip()
                    value = value.strip()

                    if value in ['True', 'False']:
                        value = value == 'True'
                    else:
                        value = float(value)

                    if key in data[current_graph]:
                        data[current_graph][key] = value

    return data



def group_by_width(parsed_data):

    # Flexible width ranges
    width_ranges = {
        "1-3": (1, 3),
        "4-6": (4, 6),
        "7-9": (7, 9),
        "10+": (10, 10000)
    }
    
    grouped_data = {key: {
        'graphs': 0, 'edges': [], 'widths': [], 'preprocess sequences heur': [],
        'total time default': [], 'total time sequences heur': [],
        'ilp time seqs heur' : [], 
        'solved default': 0, 'solved sequences heur': 0,
        'speedup_seqs': [], 'fixed vars seqs': [],
        'total time default solo': [], 'total time sequences heur solo': [],
        'graphs_common': 0,
    } for key in width_ranges}
    for graph, info in parsed_data.items():
        n,m,w = info['n'],info['m'],info['w']
        for range_label, (low, high) in width_ranges.items():
            if low <= w <= high:
                group = grouped_data[range_label]
                group['graphs'] += 1
                group['edges'].append(m)
                group['widths'].append(w)

                # Store ILP times in instances solved in all configurations
                if info['solved default'] and info['solved sequences heur']:
                    group['total time default'].append(info['total time default'])

                    group['total time sequences heur'].append(info['total time sequences heur'])
                    group['ilp time seqs heur'].append(info['ilp time seqs heur'])

                    group['graphs_common'] += 1 # Find intersection of solved instances in each safety setting

                # Preprocessing time, number of instances solved, total time for every safety setting, fixed variables
                if info['solved default']:
                    group['solved default'] += 1
                    group['total time default solo'].append(info['total time default'])
                    
                if info['solved sequences heur']:
                    group['preprocess sequences heur'].append(info['preprocess sequences heur'])
                    group['solved sequences heur'] += 1
                    group['total time sequences heur solo'].append(info['total time sequences heur'])
                    group['fixed vars seqs'].append(info['fixed vars seqs']/(w*m))

                if info['solved default'] and info['solved sequences heur']:
                    assert(info['total time sequences heur'] > 0 and info['ilp time seqs heur'] > 0)
                    #if info['total time sequences heur'] > 0 and info['ilp time seqs heur'] > 0:  # Ensure we avoid division by zero
                    speedup_seqs              = info['total time default'] / info['ilp time seqs heur']
                    speedup_seqs_with_preproc = info['total time default'] / info['total time sequences heur']
                    assert( abs(info['total time sequences heur'] - (info['ilp time seqs heur'] + info['preprocess sequences heur'])) < 0.1)
                    group['speedup_seqs'].append((speedup_seqs, speedup_seqs_with_preproc))
                
                if not info['solved default'] and info['solved sequences heur']:
                    speedup_seqs              = 300 / info['ilp time seqs heur']
                    speedup_seqs_with_preproc = 300 / info['total time sequences heur']
                    group['speedup_seqs'].append((speedup_seqs, speedup_seqs_with_preproc))
                if info['solved default'] and not info['solved sequences heur']:
                    speedup_seqs              = info['total time default'] / (300 + info['total time default']) 
                    speedup_seqs_with_preproc = info['total time default'] / (300 + info['total time default'])
                    group['speedup_seqs'].append((speedup_seqs, speedup_seqs_with_preproc))

    return grouped_data



def calculate_metrics(grouped_data):
    results = {}

    for width_range, group in grouped_data.items():
        
        results[width_range] = {
            'graphs_common': group['graphs_common'],
            'ilp_no_safety': -1,
            'ilp_safe_seqs': -1,
            'speedup_seqs': -1
        }
            
        # Calculate ILP times
        assert(len(group['total time default']) == len(group['total time sequences heur']))
        if len(group['total time default'])>0:
            results[width_range]['ilp_no_safety']  = sum(group['total time default'])  / len(group['total time default'])
            results[width_range]['ilp_safe_seqs']  = sum(group['ilp time seqs heur'])  / len(group['total time default'])

        if len(group['speedup_seqs']) > 0:
            speedup_seqs               = sum(pair[0] for pair in group['speedup_seqs']) / len(group['speedup_seqs'])
            #speedup_seqs_with_preproc = sum(pair[1] for pair in group['speedup_seqs']) / len(group['speedup_seqs'])
            results[width_range]['speedup_seqs'] = speedup_seqs

    return results



def calculate_solved(grouped_data):
    results = {}
    for width_range, group in grouped_data.items():

        results[width_range] = {
            'graphs': group['graphs'],
            'max_width': -1,
            'edges':  -1,
            'max_edges': -1,
            'preprocess_seqs': -1,
            'solved_default': -1,
            'solved_seqs': -1,
            'avg_time_default' : -1,
            'avg_time_sequences' : -1,
            'fixed_seqs': -1,
            'speedup_seqs': -1
        }

        results[width_range]['solved_default'] = group['solved default']
        results[width_range]['solved_seqs']    = group['solved sequences heur']

        # Calculate averages for preprocessing times
        if group['graphs'] > 0:
            results[width_range]['edges'] = sum(group['edges'])/group['graphs'],
            if len(group['preprocess sequences heur']) > 0:
                results[width_range]['preprocess_seqs'] = (sum(group['preprocess sequences heur']) / len(group['preprocess sequences heur']))
            results[width_range]['max_edges'] = max(group['edges'])
            results[width_range]['max_width'] = max(group['edges'])

        # Compute average total times in every safety setting
        if group['solved default'] > 0:
            results[width_range]['avg_time_default']    = sum(group['total time default solo']) / group['solved default']
        if group['solved sequences heur'] > 0:
            results[width_range]['avg_time_sequences']  = sum(group['total time sequences heur solo']) / group['solved sequences heur']

        # Calculate average of fixed vars on solved instances
        if group['solved sequences heur'] > 0:
            results[width_range]['fixed_seqs']  = 100 * sum(group['fixed vars seqs']) / group['solved sequences heur']  

        if len(group['speedup_seqs']) > 0:
            speedup_seqs               = sum(pair[0] for pair in group['speedup_seqs']) / len(group['speedup_seqs'])
            results[width_range]['speedup_seqs'] = speedup_seqs

    return results


def generate_table1(results):
    latex_code = r'''\begin{table}[]
                    \caption{Speed up metrics}
                    \begin{center}
                    \begin{tabular}{|r|r|r||r|r|r|}
                    \hline
                    & \multirow{2}{*}{width} & \multirow{2}{*}{\#graphs-all} & \multicolumn{2}{c||}{Avg.~ILP time (s)} & \multirow{2}{*}{Avg.~speedup ($\times$)} \\ \cline{4-5}
                    &                        &                           & No safety & Safe seqs.                          &                                        \\ \hline

                    \multirow{3}{*}{\rotatebox{90}{\shortstack{\textbf{Dataset}\\\textbf{name}}}}'''
 
    for width_range, metrics in results.items():
        ilp_no_safety = f"{metrics['ilp_no_safety']:.3f}" if metrics['ilp_no_safety'] != -1 else "-"
        ilp_safe_seqs = f"{metrics['ilp_safe_seqs']:.3f}" if metrics['ilp_safe_seqs'] != -1 else "-"

        speedup_seqs  = f"{metrics['speedup_seqs']:.1f}" if metrics['speedup_seqs'] != -1 else "-"

        latex_code += f"& {width_range} & {metrics['graphs_common']} & {ilp_no_safety} & {ilp_safe_seqs} & {speedup_seqs} \\\\\n"
    
    latex_code += r'''
                    \hline
                    \end{tabular}
                    \end{center}
                    \end{table}
                    '''
    return latex_code


def generate_table2(results):
    latex_code = r'''\begin{table}[]
                    \caption{Solved and fixed variables (simplified view).}
                    \begin{center}
                    \begin{tabular}{|r|r|r|r|r|r||r|r||r|}
                    \hline
                    & \multirow{2}{*}{$w$} & \multirow{2}{*}{\#g} & \multirow{2}{*}{avg $m$} & \multirow{2}{*}{preproc (s)} & \multirow{2}{*}{\shortstack{vars \\ (\%)}} & \multicolumn{2}{c|}{\#solved (avg time (s))} & \multirow{2}{*}{$\times$} \\ \cline{7-8}
                    &      &        &         &      &   & no safety & safety &  \\ \hline

                    \multirow{3}{*}{\rotatebox{90}{\shortstack{\textbf{Dataset}\\\textbf{name}}}}'''

    for width_range, metrics in results.items():
        preprocess_seqs  = f"{metrics['preprocess_seqs']:.3f}" if metrics['preprocess_seqs'] != -1 else "-"
        
        solved_default_time   = (f"{metrics['solved_default']}" if metrics['solved_default'] != -1 else "-") + " (" + (f"{metrics['avg_time_default']:.3f}"   if metrics['avg_time_default'] != -1 else "-")    + ")"
        solved_sequences_time = (f"{metrics['solved_seqs']}" if metrics['solved_seqs'] != -1 else "-")       + " (" + (f"{metrics['avg_time_sequences']:.3f}" if metrics['avg_time_sequences'] != -1 else "-")  + ")"
        
        fixed_sequences       = f"{metrics['fixed_seqs']:.1f}" if metrics['fixed_seqs'] != -1 else "-"

        speedup_seqs          = f"{metrics['speedup_seqs']:.1f}" if metrics['speedup_seqs'] != -1 else "-"

        edges_info            = ((f"{int(metrics['edges'])}") + " (" + f"{metrics['max_edges']}" + ")" ) if metrics['edges'] != -1 else "-"

        latex_code += f"& {width_range} & {metrics['graphs']} & {edges_info} & {preprocess_seqs} & {fixed_sequences} & {solved_default_time} & {solved_sequences_time} & {speedup_seqs} \\\\\n"
    
    latex_code += r'''
                    \end{tabular}
                    \end{center}
                    \end{table}
                    '''
    return latex_code



def main():

    parser = argparse.ArgumentParser(description='Process inputs.')

    parser.add_argument('-i', '--input'  , required=True, help='Input file path')

    args = parser.parse_args()

    filename    = args.input
    parsed_data = parse_input_file(filename)

    grouped_data = group_by_width(parsed_data)

    '''
    The variables `results1'' and ``latex_code1'' shoudl be ignored.
    '''
    results1      = calculate_metrics(grouped_data)
    latex_code1  = generate_table1(results1)

    results2      = calculate_solved(grouped_data)
    latex_code2   = generate_table2(results2)


    with open(filename+"1.tex", "w") as f:
        f.write(latex_code1)

    with open(filename+"2.tex", "w") as f:
        f.write(latex_code2)



if __name__ == "__main__":
    main()
