import click, sys, os, psutil, multiprocessing
from basico import *
import pandas as pd
import configparser
from basico.callbacks import create_default_handler
from multiprocessing import get_context, set_start_method
from scipy.stats import scoreatpercentile, percentileofscore
import scipy.optimize
import matplotlib.pyplot as plt
#import logging 
#logging.basicConfig(level=logging.DEBUG)

#def opener(path, DIR, flags=None):
 #   return os.open(path, flags, dir_fd=DIR)
def generalized_tanh(x, S1, S2, B):
        A = (S2 - S1) / 2
        C = (S2 + S1) / 2
        return A * np.tanh(B * x) + C

def defined_error_callback(error_obj):
        print("TimeoutError: the model creation and param estimation ran out of time. Will not continue for rest of the genes...")
        print(f'TimeoutError: {error_obj}', flush=True)

def plotPercentileWithRanks(originalPauseSums, valForNormalization, percentileval,  baseline_name, experiment_name, outputdir):
         originalPauseSums = np.sort(originalPauseSums)
         N_len = len(originalPauseSums)
         plt.clf()
         plt.plot(np.arange(1, N_len + 1), originalPauseSums, label = "Original Data", linewidth=2)
         #if logYPlotter == "Yes":
         #        plt.yscale('log')
         plt.xlabel('rank of pause sum in original data', fontsize=22)
         plt.ylabel('pause sums', fontsize=22)
         labelMax = "saturation value = %0.2f (at percentile %0.2f)" %(valForNormalization, percentileval)
         plt.plot([1,N_len+1], [valForNormalization, valForNormalization], linestyle='--', label = labelMax, linewidth=1.5)
         plt.legend()
         plt.tick_params(axis='both', which='major', labelsize=15)
         fileName = '%s_%s_percentile_%0.2f_saturationFit.pdf' %(experiment_name, baseline_name, percentileval)
         plt.savefig('%s/%s' %(outputdir, fileName), bbox_inches='tight', pad_inches=0.55)
         print('ranked paused sums with threshold at percentile %0.2f saved as %s/%s' %(percentileval, outputdir, fileName))
         # plot log Y data
         plt.clf()
         plt.plot(np.arange(1, N_len + 1), originalPauseSums, label = "Original Data", linewidth=2)
         plt.yscale('log')
         plt.xlabel('rank of pause sum in original data', fontsize=22)
         plt.ylabel('pause sums', fontsize=22)
         labelMax = "saturation value = %0.2f (at percentile %0.2f)" %(valForNormalization, percentileval)
         plt.plot([1,N_len+1], [valForNormalization, valForNormalization], linestyle='--', label = labelMax, linewidth=1.5)
         plt.legend()
         plt.tick_params(axis='both', which='major', labelsize=15)
         fileName = '%s_%s_percentile_%0.2f_saturationFit_log.pdf' %(experiment_name, baseline_name, percentileval)
         plt.savefig('%s/%s' %(outputdir, fileName), bbox_inches='tight', pad_inches=0.55)
         print('(logscale) ranked paused sums with threshold at percentile %0.2f saved as %s/%s' %(percentileval, outputdir, fileName))
         # 
def plotter_multi(originalPauseSums, percentileFitDf, baseline_name, experiment_name, outputdir, logYPlotter = "Yes"):
         N_len = len(originalPauseSums)
         #N_fitted = len(xsForSaturation)
         plt.clf()
         plt.plot(np.arange(1, N_len + 1), originalPauseSums, label = "Original Data", linewidth=2)
         for index, row in percentileFitDf.iterrows():
                 N_fitted = row['rankFitted']
                 fittedPercentile = row['percentile']
                 saturationVal = row['maxSaturation']
                 plt.plot(np.arange(1, N_fitted + 1), row['fittedVals'], label = "Fitted Data at percentile %0.1f" %(fittedPercentile), linewidth=2)
                 labelMax = "saturation (%0.1f fitted) = %0.1f (= percentile %0.1f)" %(fittedPercentile, saturationVal, percentileofscore(originalPauseSums, saturationVal))
                 plt.plot([1,len(originalPauseSums)+1], [saturationVal, saturationVal], linestyle='--', label = labelMax, linewidth=1.5)

         if logYPlotter == "Yes":
                 plt.yscale('log')
         plt.xlabel('rank of pause sum in original data', fontsize=22)
         plt.ylabel('pause sums', fontsize=22)
         plt.legend()
         plt.tick_params(axis='both', which='major', labelsize=15)
         fileName = '%s_%s_saturationFit_%s_logY.pdf' %(experiment_name, baseline_name, logYPlotter)
         plt.savefig('%s/%s' %(outputdir, fileName), bbox_inches='tight', pad_inches=0.55)
         print('saturation fits saved as %s/%s' %( outputdir, fileName))

def plotter(originalPauseSums, xsForSaturation, fittedSaturation, maxSaturation, restOfFuncEval, baseline_name, experiment_name, outputdir, logYPlotter="No"):
         N_len = len(originalPauseSums)
         N_fitted = len(xsForSaturation)
         plt.clf()
         plt.plot(np.arange(1, N_len + 1), originalPauseSums, label = "Original Data", linewidth=2)
         plt.plot(np.arange(1, N_fitted + 1), originalPauseSums[np.arange(0,N_fitted)], label = "Data For Fitting", linewidth=2)
         plt.plot(xsForSaturation, fittedSaturation, label = "Fitted Data", linewidth=2)
         if N_fitted < N_len:
                 plt.plot(np.arange(N_fitted + 1, N_len + 1), restOfFuncEval, label = "Rest of Function Evalulation", linewidth=2)
         if logYPlotter == "Yes":
                 plt.yscale('log')
         plt.xlabel('rank of pause sum in original data', fontsize=22)
         plt.ylabel('pause sums', fontsize=22)
         labelMax = "saturation (fitted curve) = %0.2f (at percentile %0.2f)" %(maxSaturation, percentileofscore(originalPauseSums, maxSaturation))
         plt.plot([1,len(originalPauseSums)+1], [maxSaturation, maxSaturation], linestyle='--', label = labelMax, linewidth=1.5)
         plt.legend()
         plt.tick_params(axis='both', which='major', labelsize=15)
         fileName = '%s_%s_percentile_%0.2f_saturationFit_%s_logY.pdf' %(experiment_name, baseline_name, (N_fitted / N_len) *100, logYPlotter)
         plt.savefig('%s/%s' %(outputdir, fileName), bbox_inches='tight', pad_inches=0.55)
         print('saturation fit for percentile %0.2f data saved as %s/%s' %((N_fitted / N_len) *100, outputdir, fileName))
# fits pause sums to hyperbolic tangent function after log-transforming them
def saturationValFromHyperTanFunc(pauseSumArr, baseline_name, experiment_name, outputdir, percentileVals = [90,95,100]):
        pauseSums = pauseSumArr[ pauseSumArr!= 0 ]
        pauseSums = np.sort(pauseSums)
        pauseSums_log10 = np.log10(pauseSums)
        N_len = len(pauseSums_log10)
        #xRanks = np.arange(1, N_len+1)
        percentileToIndices = (pd.Series(percentileVals) * 0.01 * N_len).tolist()
        listforData = []
        maxIndex = max(percentileToIndices)
        satAtMaxIndex = -1
        percentileFitData = []
        for ind in percentileToIndices:
                index = round(ind)
                xsForSaturation = np.arange(1,index + 1)
                pauseSums_log10_toFit = pauseSums_log10[np.arange(0, index)]
                pauseSums_toFit = pauseSums[np.arange(0, index)]
                try:
                        paramsSat, cv = scipy.optimize.curve_fit(generalized_tanh, xsForSaturation, pauseSums_log10_toFit, maxfev=10000)
                        S1, S2, B = paramsSat
                        A = (S2 - S1) / 2
                        C = (S2 + S1) / 2
                        maxSaturation = np.power(10,S2)
                        fittedSaturation = np.power(10, generalized_tanh(np.array(xsForSaturation), S1, S2, B))
                        restOfFuncEval = np.power(10, generalized_tanh(np.arange(index + 1, N_len + 1), S1, S2, B))
                        #plotter(pauseSums, xsForSaturation, fittedSaturation, maxSaturation, restOfFuncEval, baseline_name, experiment_name, outputdir, logYPlotter="Yes")
                        #plotter(pauseSums, xsForSaturation, fittedSaturation, maxSaturation, restOfFuncEval, baseline_name, experiment_name, outputdir, logYPlotter="No")
                        listforData.append([index, round((index/N_len)*100,1),  A, B, C, round(maxSaturation,1), round(percentileofscore(pauseSums, maxSaturation),2)])
                        percentileFitData.append([index, round((index/N_len)*100,2), maxSaturation, fittedSaturation])
                        if index == maxIndex:
                                print('will use %0.1f as saturation value corresponding to max index %d (at percentile %0.1f) provided' %(maxSaturation, maxIndex, (maxIndex/N_len)*100))
                                satAtMaxIndex = maxSaturation 
                except Exception as error:
                        print(f"Curve fit failed for index {index}: {error}")
                        listforData.append([index, round((index/N_len)*100,2), None, None, None, None, None])
                        percentileFitData.append([index, round((index/N_len)*100,2), None])
        # plot the fits with original data
        percentileFitsDf = pd.DataFrame(percentileFitData, columns = ['rankFitted', 'percentile', 'maxSaturation', 'fittedVals'])
        plotter_multi(pauseSums, percentileFitsDf, baseline_name, experiment_name, outputdir)
        plotter_multi(pauseSums, percentileFitsDf, baseline_name, experiment_name, outputdir, logYPlotter = "No")
        parameterDf = pd.DataFrame(listforData, columns = ['rankFitted','percentileRankFitted', 'A','B','C', 'saturationVal', 'percentileOfSaturationVal'])
        filename = '%s_vs_%s_fittedParameters.txt' %(experiment_name, baseline_name)
        parameterDf.to_csv('%s/%s' %(outputdir, filename), sep="\t", index=False)
        print('parameters sets written as %s/%s' %(outputdir, filename))
        return satAtMaxIndex
#remove unsuccessful fits from reports
def remove_unsuccesful_fits_from_report(filename, gene_name, baseline_name, experiment_name):
    if os.path.exists(filename):
        data_frame = pd.read_csv(filename, sep="\t")
    else:
        print("file %s not found" % filename)
        return pd.DataFrame()
    N_gene = 1 # only one gene
    new_indices = []
    count = 0
    measured_colnames = []
    fitted_colnames = []
    # first 6 * N_gene columns will be the following
    #for gene in gene_list:
    #    p_b_base_fitted = 'P_to_B_ratio_{}_Fitted'.format(baseline_name)
    #    p_b_base_measured = 'P_to_B_ratio_{}_Measured'.format(baseline_name)
    #    p_b_expr_fitted = 'P_to_B_{}_Fitted'.format(experiment_name)
    #    p_b_expr_measured = 'P_to_B_{}_Measured'.format(experiment_name)
    #    p_fc_fitted = 'P_fold_change_Fitted'
    #    p_fc_measured = 'P_fold_change_Measured'
    #indices = []
    nrow = len(data_frame)
    print("received %d parameter sets from param estimation scan" % nrow)
    for ind,row in data_frame.iterrows():
        include = True
        j = 0 # first check 4 pairs pf fitted/measured values for p_base, p_expr,  p_b_base, p_b_expr upto column index 6*N_gene
        while (j < 8 * N_gene):
            fitted_val = row[j]
            measured_val = row[j+1]
            j = j + 2
            if (abs(fitted_val - measured_val) > 0.02 * measured_val):
                include = False
                count = count + 1
                break
        if not include:
            continue
        new_indices.append(ind) # include this row as it contains valid fits
        
    #print(new_indices)
    #print(indices)
    print('>>> removed %d rows from param estimation scans' % (count))
    result_df = data_frame.loc[new_indices]
    return(result_df)


# for parsing config file containing upper and lower bound of limits
def parse_config(limitfile):
   config = configparser.ConfigParser()
   try:
        config.read(limitfile)
   except Exception as inst:
        click.echo('Error while reading the limitfile %s' %(limitfile))
        click.echo(inst)
        return {}
   section_keys = ['upper', 'lower']
   valid_section_names = ['initiation', 'pauserelease', 'elongation', 'prematuretermination', 'percentilevals']
   if len(config.sections()) == 0:
       click.echo('Error - the config file has no sections defined. Please define a config file understood by configparser (python) and refer to documentation of this package.')
       return {}
   sections_found = []
   for section in config.sections():
       if section not in valid_section_names:
           click.echo("Error - section name %s not valid. It can be either of (%s). See documentation." % (section, ",".join(valid_section_names)))
           return {}
       sections_found.append(section)
   config_dict = {}
   for section in sections_found:
        config_dict[section] = {}
        if len(config[section].keys()) == 0:
            click.echo('Error - section %s in limitfile has no fields' % (section))
            return {}
        if section == 'percentilevals':
                if 'values' not in config[section].keys():
                        click.echo('Config error: `values` not present under section `percentilevals`. Please check fomrat of `percentilevals` from vignette.')
                        return {}
                percentilevals = config.get(section,'values').split(',')
                percentilevals = [float(value.strip()) for value in percentilevals]
                if len(percentilevals) == 0 or any(value > 100 for value in percentilevals) or any(value < 1 for value in percentilevals):
                        click.echo('Config error: tried reading `values` under `percentilevals`. Either the list is empty, not separated by comma or has values outside [1-100]. Percentile values should be within [1-100] and separated by commas. Please read the vignette for syntax of the config file.')
                        return {}
                config_dict['percentilevals'] = percentilevals
        else:
                for key in config[section].keys():
                    if key not in section_keys:
                        click.echo("Error - section %s has a key %s which is not in valid keys - (%s)" %(section, key, ",".join(section_keys)))
                        return {}
                    val = config[section][key]
                    try:
                        realval = float(val)
                        config_dict[section][key] = realval
                    except Exception as inst:
                        click.echo('Error while parsing %s as a real number under section %s' %(key, section))
                        click.echo(inst)
   return config_dict

# for validating inputs
def check_inputs(basename, basefile, exprname, exprfile, limitfile, allgenesfile, norm, percentileVal, Nscan=1, algo='l'):
    click.echo('checking inputs')
    if Nscan <= 0 or Nscan > 1e7:
        click.echo('nscan cannot be zero, a negative number nor greater than 10 million. Please see help section about the usage of this parameter.')
        return False
   # if chunksz <= 0 or chunksz > 10:
   #     click.echo('chunk size (`chunksz`) cannot be zero, a negative number nor greater than 10 (not feasible). Please see help section about the usage of this parameter.')
   #     return False
    if not os.path.isfile(allgenesfile):
           click.echo('Error: the allgenesfile %s was not found' %(allgenesfile))
           return False
    if norm not in ['percentile', 'saturation']:
           click.echo("Error: unrecognized --norm %s. Can be either 'percentile' or 'saturation'. Please see `help` menu or vignette for details" %(norm))
           return False
    if percentileVal < 20 or percentileVal > 100:
            click.echo("Error: --percentileVal given as %0.2f. Must not be less than 20 or greater than 100. Please see `help` menu or vignette for details." %(percentileVal))
            return False
    algo_lower = algo.lower()
    if algo_lower not in  ['l', 'p', 'de', 'sres', 'ep', 'ga', 'gasr', 'hj', 'nl2sol', 'nm', 'prax', 'random', 'scatter', 'anneal', 'descent', 'truncated']:
            click.echo("Error: unrecognized --algo %s. Algorithm list is available in the vignette. `nl2sol` is recommended " %(algo))
            return False
    if not os.path.isfile(basefile):
        click.echo('Error: the input file for baseline was not found')
        return False
    if not os.path.isfile(exprfile):
        click.echo('Error: the input file for experiment was not found')
        return False
    if len(limitfile) !=0 and not os.path.isfile(limitfile):
        click.echo('Error: the limitfile %s was not found'  % (limitfile))
    # start reading files
    baseline_file = pd.read_csv(basefile, sep='\t')
    expr_file = pd.read_csv(exprfile, sep='\t')
    baseline_cols = baseline_file.columns.values.tolist()
    expr_cols = expr_file.columns.values.tolist()
    baseline_cols.sort()
    expr_cols.sort()
    if baseline_cols != ['body_density', 'gene', 'pause_sum']:
        click.echo("Error: corresponding columns in the baseline input file must have names 'body_density', 'gene' and 'pause_sum'. Please check the column names.")
        return False
    if expr_cols != ['body_density', 'gene', 'pause_sum']:
        click.echo("Error: corresponding columns in the experiment input file must have names 'body_density', 'gene' and 'pause_sum'. Please check the column names.")
        return False
    # check datatypes
    if baseline_file['gene'].dtype != 'O' or expr_file['gene'].dtype != 'O':
        click.echo('Error - genes in input files should be strings')
        return False
    if len(limitfile) != 0: # config file with limits have been given
        limitDict = parse_config(limitfile)
        if limitDict == {}:
            click.echo('Error - the limitfile %s had some problems. Please fix the errors.' % (limitfile))
            return False
    #force gene equivalence between two files
    baseline_genes = baseline_file['gene'].values.tolist()
    expr_genes = expr_file['gene'].values.tolist()
    baseline_genes.sort()
    expr_genes.sort()
    if baseline_genes != expr_genes:
        #print(baseline_genes)
        #print(expr_genes)
        click.echo('Error - genes in two input files dont match. Please check if files are tab separated.')
        return False
    return True

def find_rates(baseline_name, baseline_df, experiment_name, condition_df):
        rateList = []
        kpre_by_krel = 5
        for index, row in baseline_df.iterrows():
                gene = row['gene']
                normP_base = row['pause_sum']
                normB_base = row['body_density']
                exprData = condition_df.loc[condition_df['gene'] == gene]
                normP_treat = exprData["pause_sum"].values[0]
                normB_treat = exprData["body_density"].values[0]
                B_by_P_base = normB_base / normP_base
                B_by_P_treat = normB_treat / normP_treat 
                krel_FC = B_by_P_treat / B_by_P_base
                P_FC = normP_treat/normP_base # kinit FC bound 1
                kpre = B_by_P_base * 30 * kpre_by_krel ## krel * kpre/krel
                rates = [gene, normP_base, normB_base, normP_treat, normB_treat, 
                         B_by_P_base, 30*B_by_P_base, 
                         #45*B_by_P_base, 60*B_by_P_base,
                         B_by_P_treat, 30*B_by_P_treat, 
                         kpre_by_krel, kpre,
                         #45*B_by_P_treat, 60*B_by_P_treat,
                         krel_FC, P_FC, P_FC * krel_FC,
                         normB_base * 30 + kpre * normP_base,
                         normB_treat * 30 + kpre * normP_treat,
                         ]
                rateList.append(rates)

        rateDf = pd.DataFrame(rateList, columns = ['gene', 'P_normalized_%s' %(baseline_name), 'B_normalized_%s' %(baseline_name), 
                'P_normalized_%s' %(experiment_name), 'B_normalized_%s' %(experiment_name),
                'B_by_P_%s' %(baseline_name), 'pauseRelease_kelong30_%s'%(baseline_name),
                #'pauseRelease_kelong45_%s'%(baseline_name),'pauseRelease_kelong60_%s'%(baseline_name),
                'B_by_P_%s' %(experiment_name), 'pauseRelease_kelong30_%s'%(experiment_name),
                'PrematureTermination_by_PauseRelease_%s' %(baseline_name), 'PrematureTermination',
                #'pauseRelease_kelong45_%s'%(experiment_name),'pauseRelease_kelong60_%s'%(experiment_name),
                "pauseRelease_FC", "initiationRate_FC_bound_P_FC", "initiationRate_FC_bound_B_FC", 
                "initiationRate_%s_kelong30" %(baseline_name),
                #"initiationRate_%s_kelong45_termination_point1" %(baseline_name),"initiationRate_%s_kelong45_termination_1" %(baseline_name),
                "initiationRate_%s_kelong30" %(experiment_name)
                #"initiationRate_%s_kelong45_termination_point1" %(experiment_name),"initiationRate_%s_kelong45_termination_1" %(experiment_name)
                ])
        return rateDf


## NOTE: THIS FUNCTION WILL BE CALLED IN PARALLEL USING THE `MULTIPROCESSING` MODULE.   
#def execute_model(baseline_name, baseline_df, experiment_name, condition_df, index_df, Ngene, nscan, algo, thresh, kinit_limits, kelong_limits, krel_limits, kpre_limits, outputdir, geneInfo_dir, geneConstraint=False, geneConstraintObj={}):
def execute_model(baseline_name, baseline_df, experiment_name, condition_df, index_df, nscan, algo, kinit_limits, kelong_limits, krel_limits, kpre_limits, pauseMax, outputdir, geneInfo_dir):
    #click.echo('creating model for genes from index %d to %d' % (index_df, index_df + Ngene - 1))
    click.echo('creating model for gene %s' % (baseline_df.iloc[index_df]['gene']))
    # INITIALIZE MODEL
    #if geneConstraint:
    #        modelfile = "%s/compartmentModel_index_%d_to_%d_reportIndex_%d.cps" %(outputdir, index_df, index_df + Ngene - 1, geneConstraintObj["validIndex"])
    #else:
    modelfile = "%s/compartmentModel_gene_%s.cps" %(outputdir, baseline_df.iloc[index_df]['gene'])
    #if os.path.exists(modelfile): # to be removed later
    #    load_model(modelfile)
    model_notes='<body xmlns="http://www.w3.org/1999/xhtml"><h1>Compartment model</h1></body>'
    new_model(name='Compartment Model Gene_%s' % (baseline_df.iloc[index_df]['gene']), time_unit='1', quantity_unit='1', length_unit='1',area_unit='1', volume_unit='1', notes=model_notes);
    for gene_ind in range(0, 1): # just one gene
        gene_name = baseline_df.iloc[index_df + gene_ind]['gene']
        add_parameter('kelong_{}'.format(gene_name), type='global', initial_value=50)
        add_parameter('kinit_base_{}'.format(gene_name), type='global', initial_value=0.1)
        add_parameter('kinit_expr_{}'.format(gene_name), type='global', initial_value=0.1)
        add_parameter('krel_base_{}'.format(gene_name), type='global', initial_value=0.01)
        add_parameter('krel_expr_{}'.format(gene_name), type='global', initial_value=0.01)
        add_parameter('kpre_{}'.format(gene_name), type='global', initial_value=0.01) 
        # add P signal in baseline and experiment to set constraints
        #add_parameter('P_base_{}'.format(gene_name), type='assignment', expression = "Values[kinit_base_{}] / ( Values[kpre_{}] + Values[krel_base_{}] )".format(gene_name, gene_name, gene_name))
        #add_parameter('P_expr_{}'.format(gene_name), type='assignment', expression = "Values[kinit_expr_{}] / ( Values[kpre_{}] + Values[krel_expr_{}] )".format(gene_name, gene_name, gene_name))
        # add B values - B base inverse to constrain (1 / B ) <= 1 million bp
        #add_parameter('B_base_{}'.format(gene_name), type='assignment', expression = "( Values[krel_base_{}] * Values[P_base_{}] ) / Values[kelong_{}]".format(gene_name, gene_name, gene_name))
        #add_parameter('B_expr_inverse', type='assignment', expression = " Values[kelong]  / ( Values[krel_expr] * Values[P_expr] ) ")
        # add species
        add_species('P_base_{}'.format(gene_name), type='assignment', expression = "Values[kinit_base_{}] / ( Values[kpre_{}] + Values[krel_base_{}] )".format(gene_name, gene_name, gene_name))
        add_species('P_expr_{}'.format(gene_name), type='assignment', expression = "Values[kinit_expr_{}] / ( Values[kpre_{}] + Values[krel_expr_{}] )".format(gene_name, gene_name, gene_name))
        add_species(name = 'P_to_B_baseline_{}'.format(gene_name), type='assignment', expression="Values[kelong_{}] /  Values[krel_base_{}] ".format(gene_name,gene_name))
        add_species(name = 'P_to_B_expr_{}'.format(gene_name), type='assignment', expression="Values[kelong_{}] / Values[krel_expr_{}] ".format(gene_name, gene_name))
        #add_species('B_base_inverse_{}'.format(gene_name), type='assignment', expression = " Values[kelong_{}]  / ( Values[krel_base_{}] * Values[P_base_{}] ) ".format(gene_name, gene_name, gene_name))
        #add_species(name = 'P_fold_change_{}'.format(gene_name), type='assignment', expression=" ( Values[kinit_expr_{}] / Values[kinit_base_{}] ) * (  ( Values[kpre_{}] + Values[krel_base_{}] ) / ( Values[kpre_{}] + Values[krel_expr_{}] ) ) ".format(gene_name, gene_name, gene_name, gene_name, gene_name, gene_name))
    # if gene constraint - update the model to have one extra gene (index zero of baseline_df) 
    #nrow = Ngene
    gene0 =  baseline_df.iloc[index_df]['gene']
    save_model(modelfile)
    load_model(modelfile)
    click.echo('processing gene %s' % (gene0))
    # define experiments
    exp_dict = {}
    genes_with_zero_P = []
    model_index = 0
    #for index, row in baseline_df.iterrows():
    for curr_index in range(index_df, index_df + 1): # one gene only
        ### start with removing experiments
        remove_experiments()
        row = baseline_df.iloc[curr_index]
        gene = row["gene"]
        P_conc = row["pause_sum"]
        B_conc = row["body_density"]
        dictObj = {}
        dictObj['P_base'] = P_conc
        dictObj['B_base'] = B_conc
        exprData = condition_df.loc[condition_df['gene'] == gene]
        P_conc = exprData["pause_sum"].values[0]
        B_conc = exprData["body_density"].values[0]
        dictObj['P_expr'] = P_conc
        dictObj['B_expr'] = B_conc
        if dictObj['P_base'] == 0 or dictObj['P_expr'] == 0 or dictObj['B_base'] == 0 or dictObj['B_expr'] == 0:
            genes_with_zero_P.append(gene)
            click.echo("either one of P and B values for gene {} in {}/{} is zero. Skipping param estimation.".format(gene, baseline_name, experiment_name))
            continue
        #define experimental file with concentrations from P_to_B_baseline, P_to_B_expr, P_fold_change, B_fold_change from provided data
        exp_dict['[P_base_{}]'.format(gene)] = [ dictObj['P_base'] ]
        exp_dict['[P_expr_{}]'.format(gene)] = [ dictObj['P_expr'] ]
        exp_dict['[P_to_B_baseline_{}]'.format(gene)] = [ dictObj['P_base'] / dictObj['B_base'] ]
        exp_dict['[P_to_B_expr_{}]'.format(gene)] = [ dictObj['P_expr'] / dictObj['B_expr'] ]
    add_experiment('{}_{}'.format(baseline_name, experiment_name), pd.DataFrame.from_dict(exp_dict, orient='columns'), file_name='{}/inputFile_{}_{}_gene_{}.txt'.format(geneInfo_dir, baseline_name, experiment_name, baseline_df.iloc[index_df]['gene'] ))
    # start param fitting
    param_fit_items = []
    fitConstraints = []
    #for gene in baseline_df['gene']:
    for curr_index in range(index_df, index_df + 1): # one gene only
        # TBD: add constraint - both P_baseline and P_expr <= 2
        #print('fitting params for gene %s' %(gene))
        gene_name = baseline_df.iloc[curr_index]['gene']
        #kinit
        dictObj = {}
        dictObj['name'] = 'Values[kinit_base_{}]'.format(gene_name)
        #dictObj['name'] = 'kinit_{}'.format(gene_name)
        dictObj['lower'] = kinit_limits[0]
        dictObj['upper'] = kinit_limits[1]
        dictObj['start'] = kinit_limits[0]
        #dictObj['affected'] = []
        param_fit_items.append(dictObj)
        dictObj = {}
        dictObj['name'] = 'Values[kinit_expr_{}]'.format(gene_name)
        #dictObj['name'] = 'kinit_{}'.format(gene_name)
        dictObj['lower'] = kinit_limits[0]
        dictObj['upper'] = kinit_limits[1]
        dictObj['start'] = kinit_limits[0]
        #dictObj['affected'] = ['{}_{}'.format(experiment_name,gene_name)]
        param_fit_items.append(dictObj)
        #krel
        dictObj = {}
        dictObj['name'] = 'Values[krel_base_{}]'.format(gene_name)
        #dictObj['name'] = 'krel_{}'.format(gene_name)
        dictObj['lower'] = krel_limits[0]
        dictObj['upper'] = krel_limits[1]
        dictObj['start'] = krel_limits[0]
        #dictObj['affected'] = []
        param_fit_items.append(dictObj)
        dictObj = {}
        dictObj['name'] = 'Values[krel_expr_{}]'.format(gene_name)
        #dictObj['name'] = 'krel_{}'.format(gene_name)
        dictObj['lower'] = krel_limits[0]
        dictObj['upper'] = krel_limits[1]
        dictObj['start'] = krel_limits[0]
        #dictObj['affected'] = ['{}_{}'.format(experiment_name, gene_name)]
        param_fit_items.append(dictObj)
        # kpre
        dictObj = {}
        dictObj['name'] = 'Values[kpre_{}]'.format(gene_name)
        #dictObj['name'] = 'kpre_{}'.format(gene_name)
        dictObj['lower'] = kpre_limits[0]
        dictObj['upper'] = kpre_limits[1]
        dictObj['start'] = kpre_limits[0]
        #dictObj['affected'] = []
        param_fit_items.append(dictObj)
        #kelong
        dictObj = {}
        dictObj['name'] = 'Values[kelong_{}]'.format(gene_name)
        #dictObj['name'] = 'kelong_{}'.format(gene_name)
        dictObj['lower'] = kelong_limits[0]
        dictObj['upper'] = kelong_limits[1]
        dictObj['start'] = kelong_limits[0]
        #dictObj['affected'] = []
        #print(dictObj)
        param_fit_items.append(dictObj)
        # deal with fit constraints
        #dictObj = {}
        #dictObj['name'] = 'Values[P_base_{}]'.format(gene_name)
        #dictObj['lower'] = 0
        #dictObj['upper'] = pauseMax
        #fitConstraints.append(dictObj)
        #dictObj = {}
        #dictObj['name'] = 'Values[P_expr_{}]'.format(gene_name)
        #dictObj['lower'] = 0
        #dictObj['upper'] = pauseMax
        #fitConstraints.append(dictObj)
       # dictObj = {}
       # dictObj['name'] = 'Values[B_base_inverse_{}]'.format(gene_name)
       # dictObj['lower'] = 0
       # dictObj['upper'] = 1e6
        fitConstraints.append(dictObj)
    #print(param_fit_items)
    #print(param_fit_items)
    set_fit_parameters(param_fit_items,model=get_current_model())
    set_fit_constraints(fitConstraints)
    algo = algo.lower()
    create_default_handler(max_time=20)
    if (algo == "l"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Levenberg - Marquardt', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Levenberg - Marquardt', randomize_start_values=True, update_model=True)
    elif (algo == 'p'):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Particle Swarm', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Particle Swarm', randomize_start_values=True, update_model=True)
    elif (algo == "de"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Differential Evolution', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Differential Evolution', randomize_start_values=True, update_model=True) 
    elif (algo == "sres"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Evolution Strategy (SRES)', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Evolution Strategy (SRES)', randomize_start_values=True, update_model=True)
    elif (algo == "ep"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Evolutionary Programming', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Evolutionary Programming', randomize_start_values=True, update_model=True)
    elif (algo == "ga"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Genetic Algorithm', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Genetic Algorithm', randomize_start_values=True, update_model=True)
    elif (algo == "gasr"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Genetic Algorithm SR', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Genetic Algorithm SR', randomize_start_values=True, update_model=True)  
    elif (algo == "hj"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Hooke & Jeeves', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Hooke & Jeeves', randomize_start_values=True, update_model=True)   
    elif (algo == "nl2sol"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'NL2SOL', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='NL2SOL', randomize_start_values=True, update_model=True)   
    elif (algo == "nm"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Nelder - Mead', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Nelder - Mead', randomize_start_values=True, update_model=True)   
    elif (algo == "prax"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Praxis', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Praxis', randomize_start_values=True, update_model=True)   
    elif (algo == "random"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Random Search', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Random Search', randomize_start_values=True, update_model=True)   
    elif (algo == "scatter"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Scatter Search', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Scatter Search', randomize_start_values=True, update_model=True)
    elif (algo == "anneal"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Simulated Annealing', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Simulated Annealing', randomize_start_values=True, update_model=True)   
    elif (algo == "descent"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Steepest Descent', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Steepest Descent', randomize_start_values=True, update_model=True)
    elif (algo == "truncated"):
            set_task_settings(T.PARAMETER_ESTIMATION, settings={'method':'Truncated Newton', 'randomize_start_values': True, 'scheduled': True})
            run_parameter_estimation(method='Truncated Newton', randomize_start_values=True, update_model=True)  
    else:
            raise Exception("unrecognized --algo flag %s provided. Please see help section or vignette." % (algo))
    #save_model(modelfile)
    #load_model(modelfile)
    # >>>> 'repeat' task for scan with Nscan steps
    numScanSteps = nscan
    bodyDef = []
    reportColumnNames = []
    ind = 0
    for curr_index in range(index_df, index_df + 1): # only one gene
        gene = baseline_df.iloc[curr_index]['gene']
        P_base_fitted = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Fitted Value".format(baseline_name, experiment_name, ind)
        colName = "P_{}_Fitted".format(baseline_name)
        bodyDef.append(P_base_fitted)
        reportColumnNames.append(colName)
        P_base_measured = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Measured Value".format(baseline_name, experiment_name, ind)
        colName = "P_{}_Measured".format(baseline_name)
        bodyDef.append(P_base_measured)
        reportColumnNames.append(colName)
        ind = ind + 1 # increase the index
        
        P_expr_fitted = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Fitted Value".format(baseline_name, experiment_name, ind)
        colName = "P_{}_Fitted".format(experiment_name)
        bodyDef.append(P_expr_fitted)
        reportColumnNames.append(colName)
        P_expr_measured = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Measured Value".format(baseline_name, experiment_name, ind)
        colName = "P_{}_Measured".format(experiment_name)
        bodyDef.append(P_expr_measured)
        reportColumnNames.append(colName)
        ind = ind + 1 # increase the index

        P_to_B_baseline_fitted = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Fitted Value".format(baseline_name, experiment_name, ind)
        colName = "P_to_B_ratio_{}_Fitted".format(baseline_name)
        bodyDef.append(P_to_B_baseline_fitted)
        reportColumnNames.append(colName)

        P_to_B_baseline_measured = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Measured Value".format(baseline_name, experiment_name, ind)
        colName = "P_to_B_ratio_{}_Measured".format(baseline_name)
        bodyDef.append(P_to_B_baseline_measured)
        reportColumnNames.append(colName)
        ind = ind + 1 # increase the index
        
        P_to_B_expr_fitted = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Fitted Value".format(baseline_name, experiment_name, ind)
        colName = "P_to_B_ratio_{}_Fitted".format(experiment_name)
        bodyDef.append(P_to_B_expr_fitted)
        reportColumnNames.append(colName)
        P_to_B_expr_measured = "CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,ParameterGroup=Experiment Set,ParameterGroup={}_{},Vector=Fitted Points[{}],Reference=Measured Value".format(baseline_name, experiment_name, ind)
        colName = "P_to_B_ratio_{}_Measured".format(experiment_name)
        bodyDef.append(P_to_B_expr_measured)
        reportColumnNames.append(colName)
        ind = ind + 1 # increase the index

    for curr_index in range(index_df, index_df + 1): # only one gene
        gene = baseline_df.iloc[curr_index]['gene']
        kelong =  'Values[kelong_{}]'.format(gene)
        colName = "elongationRate"
        bodyDef.append(kelong)
        reportColumnNames.append(colName)
        
        kinit_base =  'Values[kinit_base_{}]'.format(gene)
        colName = "initiationRate_{}".format(baseline_name)
        bodyDef.append(kinit_base)
        reportColumnNames.append(colName)

        kinit_expr =  'Values[kinit_expr_{}]'.format(gene)
        colName = "initiationRate_{}".format(experiment_name)
        bodyDef.append(kinit_expr)
        reportColumnNames.append(colName)

        kpre =  'Values[kpre_{}]'.format(gene)
        colName = "prematureTerminationRate"
        bodyDef.append(kpre)
        reportColumnNames.append(colName)

        krel_base =  'Values[krel_base_{}]'.format(gene)
        colName = "pauseReleaseRate_{}".format(baseline_name)
        bodyDef.append(krel_base)
        reportColumnNames.append(colName)

        krel_expr =  'Values[krel_expr_{}]'.format(gene)
        colName = "pauseReleaseRate_{}".format(experiment_name)
        bodyDef.append(krel_expr)
        reportColumnNames.append(colName)

    bestObjVal = 'CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Value'
    bodyDef.append(bestObjVal)
    reportColumnNames.append('objective_func_bestval')
    # set scan items
    #print("estimating params")
    #    # define output report and assign it to a textfile
    reportName = 'ParamEstimation_report'
    #remove_report(reportName)
    add_report(reportName, task=T.STEADY_STATE, is_table=False, header=reportColumnNames, body=bodyDef)
    repeatedParamDir = '%s/repeatedParamScans' %(os.path.abspath(outputdir))
    reportFile = '%s/repeatedParamEstimates_gene_%s.txt' % (repeatedParamDir, baseline_df.iloc[index_df]['gene'])
    assign_report(reportName, task = T.SCAN, filename = reportFile, append=False)
    set_scan_items([{'type':'repeat', 'num_steps': numScanSteps, 'output during subtask execution': False}], subtask = T.PARAMETER_ESTIMATION)
    set_scan_settings(settings={'update_model': True, 'scheduled': True,
            'subtask': T.PARAMETER_ESTIMATION,
            'output_during_subtask': False,
            'continue_from_current_state': False,
            'continue_on_error': True, 'scan_items': get_scan_items()}) 
    save_model(modelfile) # to preserve the model state 
    load_model(modelfile) # try to reproduce copasi GUI load 
    #run_scan(settings={'update_model': True, 'scheduled': True, 
    #        'subtask': T.PARAMETER_ESTIMATION, 
    #        'output_during_subtask': False, 
    #        'continue_from_current_state': False,
    #        'continue_on_error': True, 'scan_items': get_scan_items()})
    run_scan(settings = get_scan_settings())
    # add logic to filter the unsuccessful fits
    save_model(modelfile) # to preserve the model state 
    final_df = remove_unsuccesful_fits_from_report(reportFile, baseline_df.iloc[index_df]['gene'], baseline_name, experiment_name)
    reportOutFile = '{}/repeatedParamEstimates_gene_{}_onlyvalidfits.txt'.format(repeatedParamDir, baseline_df.iloc[index_df]['gene'])
    if not final_df.empty:
            #valid_gene_list.append(gene)
            final_df.to_csv(reportOutFile, sep="\t", index=False)
            click.echo(">>> saved valid fits to {}".format(reportOutFile))
            return True
    else:
            click.echo("No valid fits were found for gene %s. Please increase `nscan` (--nscan) if you are worried. See vignette for more details." % (baseline_df.iloc[index_df]['gene']))
            return False # NA value

# this function will parse a report file and get all the parts
def get_linked_params(validReportFile, reportInd, chunksz, indexStart, indexEnd, baseline_name, experiment_name):
        if not os.path.isfile(validReportFile):
                click.echo("file {} not found".format(validReportFile))
                return pd.DataFrame()
        #df = pd.DataFrame([], columns = ['gene', 'rateName', 'value', 'condition', 'paramSet', 'reportIndex'])
        parameterSets = pd.read_csv(validReportFile, sep="\t")
        columnList = parameterSets.columns.values.tolist()
        #firstRateIndices = [i for i in range(0, len(columnList), chunksz)]
        firstRates = columnList[0:(chunksz*6):6] # take every gene's elongation rate after six steps - as there are six rates
        geneList = [arr.split('_')[1] for arr in firstRates]
        listOfConnectedRates = []
        #print(geneList)
        for ind, row in parameterSets.iterrows():
                for gene in geneList:
                        kelong = row["elongationRate_{}".format(gene)]
                        kpre = row["prematureTerminationRate_{}".format(gene)]
                        kinit_base = row["initiationRate_{}_{}".format(gene,baseline_name)]
                        kinit_expr = row["initiationRate_{}_{}".format(gene,experiment_name)]
                        krel_base = row["pauseReleaseRate_{}_{}".format(gene,baseline_name)]
                        krel_expr = row["pauseReleaseRate_{}_{}".format(gene,experiment_name)]
                        listOfConnectedRates.append([gene, kinit_base, kinit_expr, krel_base, krel_expr,  kelong, kpre, ind, reportInd])
                       # listOfConnectedRates.append([gene, "elongation", kelong, baseline_name, ind, reportInd])
                       # listOfConnectedRates.append([gene, "prematureTermination", kpre, baseline_name, ind, reportInd])
                       # listOfConnectedRates.append([gene, "initiation", kinit_base, baseline_name, ind, reportInd])
                       # listOfConnectedRates.append([gene, "initiation", kinit_expr, experiment_name, ind, reportInd])
                       # listOfConnectedRates.append([gene, "pauseRelease", krel_base, baseline_name, ind, reportInd])
                       # listOfConnectedRates.append([gene, "pauseRelease", krel_expr, experiment_name, ind, reportInd])

        #df = pd.DataFrame(listOfConnectedRates, columns = ['gene', 'rateName', 'value', 'condition', 'paramSet', 'reportIndex'])
        df = pd.DataFrame(listOfConnectedRates, columns = ['gene', 'initiation_{}'.format(baseline_name), 'initiation_{}'.format(experiment_name), 'pauseRelease_{}'.format(baseline_name), 'pauseRelease_{}'.format(experiment_name), "elongation","prematureTermination", 'paramSet', 'reportIndex'])
        return df

@click.group() # ran before any subcommand is executed
def cli():
    pass

# for param estimation
@cli.command()
@click.option('--basename', default='baseline', help='the identifier of the baseline condition. Will be used as suffix to relevant column names in output files', show_default=True)
@click.option('--basefile', required=True, help="path of file containing baseline body and pause densities")
#in bed6 format containing genes and densities in their body and pause regions. The 'score' column (column 5) contains the si")
@click.option('--exprname', default='experiment', help='the identifier for the experimental condition. Will be used as a suffix to relevant column names in output files', show_default=True)
@click.option('--exprfile', required=True, help="path of file containing baseline body and pause densities")
@click.option('--outputdir', required=True, help="directory to which the results will be written. For instance, paramFitDir_<experiment_name>. Say, 'paramFitDir_ZNF143_degron'.")
#@click.option('--onlyvalidfit', required=False, default=True, help="only save parameters for genes having valid fit", show_default=True)
@click.option('--conf', required=False, default="", help="only use this option if you are interested in changing the default parameter settings. This is path to a config file understood by ConfigParser() in configparser package containing the percentile vals or lower and upper limits of parameters to fit. Note - empty config files are not accepted. Valid section names are 'initiation', 'pauserelease', 'elongation', 'prematuretermination', 'percentilevals'. An example file can be found in this vignette: http://guertinlab.cam.uchc.edu/compartMentModel_vignette_RM/compartmentModel_vignette.html.", show_default=True)
#@click.option('--nscan', type=int, required=False, default=500, help="a number N denoting the number of param estimation runs to call for every gene in the input file. If not set, parameter estimation will be called ten times for every gene.", show_default=True)
#@click.option('--algo', default='l', required=False, help="algorithm for parameter estimation. Options - 'L' (Levenberg-Marquardt; default) and 'nl2sol' (NL2SOL). See vignette for other algorithms available.", show_default=True)
@click.option('--allgenesfile', help='file containing pause sums of relevant genes to be considered for calculating normalization value. The file should be tab-separated, containing atleast two columns for `gene` and `pause_sum`. See vignette for more details about how this file is defined on body density var/mean criteria.',required=True)
@click.option('--norm', default='percentile', required=False, help="normalization method to scale data to occupancy values. Can be 'percentile' OR 'saturation'. If 'percentile', the `percentileVal` pause sum from allgenesfile data (`allgenesfile`) will be used. If 'saturation', all the pause sums in `allgenesfile` data are fitted to hyperbolic tangent function and the saturation value from the fit is used for normalization. See vignette for details.", show_default=True)
@click.option('--percentileval', type=float, default = 90, help="if `norm` is 'percentile', then this percentile of pause sum from allgenes data (`allgenesfile`) will be used for normalization", required=False, show_default=True)
#@click.option('--overwrite/--no-overwrite', default=False, help='flag to overwrite a model file created previously', show_default=True)
#@click.option('--thresh', type=float, required=False, default=5.0, help="threshold percentage of genes that don't have valid parameter fits for global ratios of their body densities. Parameter sets having invalid fits for more genes than this threshold are not included in the final output. See vignette for more details.", show_default=True)
#@click.option('--chunksz', type=int, required=False, default=5, help="divide gene set into chunks of this size. The model will be called on each chunk separately. See vignette for more details.", show_default=True)
#def estparam(basename, basefile, exprname, exprfile, outputdir, conf, nscan, algo, allgenesfile, norm, percentileval):
def estparam(basename, basefile, exprname, exprfile, outputdir, conf, allgenesfile, norm, percentileval):
    ''' Estimate params kinit, krel, kelong, kpre of compartment model. See Dutta et al Genome Res, 2023 for more details.
    '''
    if not check_inputs(basename, basefile, exprname, exprfile, conf, allgenesfile, norm, percentileval, Nscan=2, algo='l'):
    #if not check_inputs(basename, basefile, exprname, exprfile, conf, allgenesfile, norm, percentileval):
            sys.exit('please fix the error with inputs and try again.')
    else:
            baseline_df = pd.read_csv(basefile, sep='\t')
            condition_df = pd.read_csv(exprfile, sep='\t')
            allgenes_df = pd.read_csv(allgenesfile, sep="\t")
            #create output directory
            if not os.path.isdir(outputdir):
                click.echo('creating output directory %s' %(outputdir))
                os.makedirs(outputdir)

            ## normalization and saturation logic
            percentileVals = [90,95,100]
            if len(conf) != 0: # config file with limits have been given
                limitDict = parse_config(conf)
                if 'percentilevals' in limitDict.keys():
                        percentileVals = limitDict['percentilevals']
            if 'pause_sum' not in allgenes_df.columns:
                    sys.exit('Error: column `pause_sum` doesnt exist in %s file. Please fix error with inputs and try again. See help and vignette for details' %(allgenesfile))
            if len(allgenes_df) < len(baseline_df):
                    sys.exit('Error: number of rows in `allgenesfile` is less than number of rows in input files. Please check provided gene sets.')
            pForNormalization = -1
            if norm == 'percentile':
                    allgenes_df =  allgenes_df.loc[allgenes_df['pause_sum'] > 0]
                    pForNormalization = scoreatpercentile(allgenes_df['pause_sum'], percentileval)
                    plotPercentileWithRanks(allgenes_df['pause_sum'], pForNormalization, percentileval,  basename, exprname, outputdir)
                    click.echo('pause sum of %0.2f (at percentile %0.2f) will be used to normalize data' %(pForNormalization, percentileval))
            else: # norm == 'saturation'
                     #allgenes_df =  allgenes_df.loc[allgenes_df['pause_sum'] > 0]
                     pForNormalization = saturationValFromHyperTanFunc(allgenes_df['pause_sum'].values, basename, exprname, outputdir, percentileVals)
                     if pForNormalization == -1:
                             sys.exit('saturation curve fit failed. Please check if provided data has reasonable number of genes and densities are computed properly. You can alternatively use `--norm=percentile` for calculating a pause sum for normalization. Please see --help menu or vignette.')
                     #click.echo('From curve fit, saturation value of %0.2f will be used to normalize data' %(pForNormalization))
            baseline_df['pause_sum'] = baseline_df['pause_sum'] / pForNormalization
            baseline_df['body_density'] = baseline_df['body_density'] / pForNormalization
            condition_df['pause_sum'] = condition_df['pause_sum'] / pForNormalization
            condition_df['body_density'] = condition_df['body_density'] / pForNormalization
            # rest of the logic
            baseline_name = basename
            experiment_name = exprname
            # discard zero P and B
            copy_baseline = baseline_df
            copy_condition = condition_df
            baseline_df = baseline_df.loc[(copy_baseline['pause_sum'] > 0) & (copy_baseline['body_density'] > 0) & (copy_condition['pause_sum'] > 0) & (copy_condition['body_density'] > 0)]
            condition_df = condition_df.loc[(copy_baseline['pause_sum'] > 0) & (copy_baseline['body_density'] > 0) & (copy_condition['pause_sum'] > 0) & (copy_condition['body_density'] > 0)]
            baseline_df = baseline_df.reset_index(drop=True)
            condition_df = condition_df.reset_index(drop=True)
            removedGenes = set(copy_baseline['gene']).difference(baseline_df['gene'])
            if len(removedGenes) != 0:
                click.echo(">> removed %d genes from run as they have zero pause sum or body density" % (len(removedGenes)))
            genes_with_zero_P = list(removedGenes)
            df_genes_zeroP = pd.DataFrame(genes_with_zero_P, columns=['Values'])
            df_genes_zeroP.to_csv("{}/GeneswithZeroPauseSums_or_BodyAvg.txt".format(outputdir), index=False, header=False)
            click.echo('>> Genes with pause sum or body avg as zero in either condition have been written to {}'.format("{}/GeneswithZeroPauseSums_or_BodyAvg.txt".format(outputdir)))
            
           # dir_descrp = os.open(outputdir, os.O_RDONLY)
           # if len(genes_with_zero_P) != 0:
           #         zeroPauseFile = "{}/GeneswithZeroPauseSums_or_BodyAvg.txt".format(outputdir)
           #         with open(zeroPauseFile, 'w') as f:
           #                 for gene in genes_with_zero_P:
           #                         print("%s" % gene, file = f)
           #         click.echo('>> Genes with pause sum or body avg as zero in either condition have been written to {}'.format(zeroPauseFile))
           # os.close(dir_descrp)
            ## discard genes with Pause Sum > 1
            copy_baseline = baseline_df
            copy_condition = condition_df
            baseline_df = baseline_df.loc[(copy_baseline['pause_sum'] <= 1) & (copy_condition['pause_sum'] <= 1) ]
            condition_df = condition_df.loc[(copy_baseline['pause_sum'] <= 1) & (copy_condition['pause_sum'] <= 1)]
            baseline_df = baseline_df.reset_index(drop=True)
            condition_df = condition_df.reset_index(drop=True)
            removedGenes = set(copy_baseline['gene']).difference(baseline_df['gene'])
            if len(removedGenes) != 0:
                click.echo(">> removed %d genes from run as their pause sum exceeds the saturation value being used" % (len(removedGenes)))
            genes_with_P_larger_1 = list(removedGenes)
            df_genes_with_P_larger_1 = pd.DataFrame(genes_with_P_larger_1,  columns=['Values'])
            df_genes_with_P_larger_1.to_csv("{}/GeneswithPauseSums_ExceedSaturationValue.txt".format(outputdir), index=False, header=False)
            click.echo('>> Genes with pause sum exceeding saturation values, written to {}'.format("{}/GeneswithPauseSums_ExceedSaturationValue.txt".format(outputdir)))
            commandUsed = " ".join(sys.argv) 
            saveCmdFile = "{}/commandUsed.txt".format(outputdir)
            print(">> Received the command:\n{}\n>>".format(commandUsed))
            print("writing command information to {}".format(saveCmdFile))
            with open(saveCmdFile, 'w') as f:
                print(commandUsed, file=f)
            geneInfo_dir = '%s/geneInfo' %(outputdir)
           # if not os.path.isdir(geneInfo_dir):
           #     click.echo('creating dir %s for storing empirical ratios for each gene' %(geneInfo_dir))
           #     os.makedirs('%s' %(geneInfo_dir)) 
           # N = len(baseline_df) 
           # multiples = int(N / chunksz)
           # remainders = int(N % chunksz)
            # set limits for param estimation
            kinit_limits = [0.01,1]
            kelong_limits = [30,60]
            krel_limits = [0.01,1]
            kpre_limits = [0.01,1]
            pauseMax = 1
            #if len(conf) != 0: # config file with limits have been given
            #    limitDict = parse_config(conf)
            #    if 'initiation' in limitDict.keys():
            #        if 'lower' in limitDict['initiation'].keys():
            #            kinit_limits[0] = limitDict['initiation']['lower']
            #        if 'upper' in limitDict['initiation'].keys():
            #            kinit_limits[1] = limitDict['initiation']['upper']
            #        if kinit_limits[0] > kinit_limits[1]:
            #            click.echo('Error - the lower limit for initiation - %f is larger than the upper limit %f' % (kinit_limits[0], kinit_limits[1]))
            #            sys.exit('please fix the error with limits file and try again.')
            #    if 'pauserelease' in limitDict.keys():
            #        if 'lower' in limitDict['pauserelease'].keys():
            #            krel_limits[0] = limitDict['pauserelease']['lower']
            #        if 'upper' in limitDict['pauserelease'].keys():
            #            krel_limits[1] = limitDict['pauserelease']['upper']
            #        if krel_limits[0] > krel_limits[1]:
            #            click.echo('Error - the lower limit for pauserelease - %f is larger than the upper limit %f' % (krel_limits[0], krel_limits[1]))
            #            sys.exit('please fix the error with limits file and try again.')
            #    if 'elongation' in limitDict.keys():
            #        if 'lower' in limitDict['elongation'].keys():
            #            kelong_limits[0] = limitDict['elongation']['lower']
            #        if 'upper' in limitDict['elongation'].keys():
            #            kelong_limits[1] = limitDict['elongation']['upper']
            #        if kelong_limits[0] > kelong_limits[1]:
            #            click.echo('Error - the lower limit for elongation - %f is larger than the upper limit %f' % (kelong_limits[0], kelong_limits[1]))
            #            sys.exit('please fix the error with limits file and try again.')
            #    if 'prematuretermination' in limitDict.keys():
            #        if 'lower' in limitDict['prematuretermination'].keys():
            #            kpre_limits[0] = limitDict['prematuretermination']['lower']
            #        if 'upper' in limitDict['prematuretermination'].keys():
            #            kpre_limits[1] = limitDict['prematuretermination']['upper']
            #        if kpre_limits[0] > kpre_limits[1]:
            #            click.echo('Error - the lower limit for prematuretermination - %f is larger than the upper limit %f' % (kpre_limits[0], kpre_limits[1]))
            #            sys.exit('please fix the error with limits file and try again.')
                #if 'maxpause' in limitDict.keys():
                #    if 'value' in limitDict['maxpause'].keys():
                #        val = float(limitDict['maxpause']['value'])
                #        if val > 3 or val <= 0:
                #            click.echo('Error - max Pause value cannot be non-positive number or greater than 3. Value under `maxpause` in config file was given as %f' %(val))
                #            sys.exit('please fix the error with limits file and try again.')
                #        pauseMax = val
            
            
           # repeatedParamDir = '%s/repeatedParamScans' %(outputdir)
           # if not os.path.isdir(repeatedParamDir):
           #     click.echo('creating directory %s' %(repeatedParamDir))
           #     os.makedirs(repeatedParamDir)
            rateDataDf = find_rates(baseline_name, baseline_df, experiment_name, condition_df)
            outfilename = "{}/paramEstimationResults_{}_vs_{}.txt".format(outputdir, experiment_name, baseline_name)
            rateDataDf.to_csv(outfilename, sep="\t", index=False)
            print("> Estimated rates saved as {}".format(outfilename))

 
            ## OLD method which invoked COPASI to execute the param estimation
            ## report on the algo used
            #algo = algo.lower()
            #if (algo == "l"):
            #        click.echo("param estimation will be run using Levenberg-Marqaurdt algorithm.")
            #elif (algo == "p"):
            #        click.echo("param estimation will be run using Particle-Swarm algorithm.")
            #elif (algo == "de"):
            #        click.echo("param estimation will be run using Differential Evolution algorithm.")
            #elif (algo == "sres"):
            #        click.echo("param estimation will be run using Evolution Strategy (SRES) algorithm.")
            #elif (algo == "ep"):
            #        click.echo("param estimation will be run using Evolutionary Programming method.")
            #elif (algo == "ga"):
            #        click.echo("param estimation will be run using Genetic Algorithm method.")
            #elif (algo == "gasr"):
            #        click.echo("param estimation will be run using Genetic Algorithm SR method.")
            #elif (algo == "hj"):
            #        click.echo("param estimation will be run using Hooke & Jeeves method.")
            #elif (algo == "nl2sol"):
            #        click.echo("param estimation will be run using NL2SOL method.")
            #elif (algo == "nm"):
            #        click.echo("param estimation will be run using Nelder - Mead method.")
            #elif (algo == "prax"):
            #        click.echo("param estimation will be run using Praxis method.")
            #elif (algo == "random"):
            #        click.echo("param estimation will be run using Random Search method.")
            #elif (algo == "scatter"):
            #        click.echo("param estimation will be run using Scatter Search method.")
            #elif (algo == "anneal"):
            #        click.echo("param estimation will be run using Simulated Annealing method.")
            #elif (algo == "descent"):
            #        click.echo("param estimation will be run using Steepest Descent method.")
            #elif (algo == "truncated"):
            #        click.echo("param estimation will be run using Truncated Newton method.")    
            #else:
            #        raise sys.exit("Error - unrecognized --algo flag %s detected. Please see help section" % (algo))

            ## start with setting cores and num of processes to be spawned
            #ncore = 0
            #if 'SLURM_CPUS_PER_TASK' in os.environ:
            #    ncore = int(os.environ['SLURM_CPUS_PER_TASK'])
            #else:
            #    ncore = psutil.cpu_count(logical=False)
            ##click.echo('>> received Nscan as %d. Chunksize is %d. Deduced Ncore as %d.' % (nscan, chunksz, ncore))
            #click.echo('>> received Nscan as %d. Deduced Ncore as %d.' % (nscan, ncore))
            ## functionSignature: execute_model(baseline_name, baseline_df, experiment_name, condition_df, index_df, Ngene, nscan, thresh, kinit_limits, kelong_limits, krel_limits, kpre_limits, outputdir)
            #
            ### to prevent python from stalling, referenced from: https://pythonspeed.com/articles/python-multiprocessing/
            ##set_start_method("spawn") 
            ##with get_context("spawn").Pool(ncore) as pool:
            #with multiprocessing.Pool(ncore) as pool:
            #    results = pool.starmap(execute_model, [(baseline_name, baseline_df, experiment_name, condition_df, index_df, nscan, algo, kinit_limits, kelong_limits, krel_limits, kpre_limits, pauseMax, outputdir, geneInfo_dir) for index_df, _ in baseline_df.iterrows()]) #error_callback=defined_error_callback)
            #    #results.get(timeout=(2.4*len(baseline_df)))
            #click.echo("Param estimation finished. Collecting valid fits...")
            #Nanchor = 1
            #noValidFits = []
            #noValidCount = 0
            ### parse report files generated and create a manageable file
            #df_list = []
            #outFile = "{}/combinedResults.txt".format(outputdir)
            #for geneName in baseline_df['gene']:
            #    #geneName = row["gene"]
            #    reportOutFile = '{}/repeatedParamEstimates_gene_{}_onlyvalidfits.txt'.format(repeatedParamDir, geneName)
            #    if os.path.exists(reportOutFile):
            #        reportDf = pd.read_csv(reportOutFile, sep="\t")
            #        reportDf = reportDf.assign(gene = geneName)
            #        df_list.append(reportDf)
            #    else:
            #        noValidFits.append(geneName)
            ## write the combined rates
            #if len(df_list) != 0:
            #    combined_df = pd.concat(df_list)
            #    outfilename = "{}/paramEstimationResults_{}_vs_{}_combined.txt".format(outputdir, experiment_name, baseline_name)
            #    combined_df.to_csv(outfilename, sep="\t", index=False)
            #    print("> Valid fits are saved as {} ".format(outfilename))
            #else:
            #        click.echo("No valid fits found. Please increase nscan or change upper/lower limits on parameters.")
            ## save genes with no valid fits
            #if len(noValidFits) != 0:
            #        #collapseList = sum(noValidFits, []) # flatten to a single column
            #        outObj = pd.DataFrame.from_dict({'genes_no_valid_fit': noValidFits})
            #        outFile = "{}/NoValidFits.txt".format(outputdir)
            #        click.echo(">>> {} genes with no valid fits will be written to {}".format(len(noValidFits), outFile))
            #        outObj.to_csv(outFile, sep="\t", index=False, header=True)

