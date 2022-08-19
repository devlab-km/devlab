import pandas as pd
import json
import os
from plotnine import *
from statsmodels.stats.proportion import proportion_confint

class covidSeqReport:
    '''Class to handle the data collection from a single sample's combined report
    from the Qiagen Software output'''

    def __init__(self, report_filename, worksheet_filename):
        try:
            self.data = self.open_json(report_filename)
        except:
            print('Missing Report')
        try:
            self.qc = self.data['data']['qc_for_sequencing_reads']
        except:
            print('Missing Sequencing Reads QC')
        try:
            self.summary = self.data['data']['reads_summary']
        except:
            print('Missing Reads Summary')
        try:
            self.trim = self.data['data']['trim_reads']
        except:
            print('Missing Trim Reads')
        try:
            self.map_summary = self.data['data']['read_mapping_summary']
        except:
            print('Missing SCOV2 Mapping Summary')
        try:
            self.duplicates = self.data['data']['duplicated_mapped_reads']
        except:
            print('Missing Duplicates Report')
        try:
            self.qc_targeted = self.data['data']['qc_for_targeted_sequencing']
        except:
            print('Missing Targeted QC Report')
        try:
            self.sample_summary = self.open_worksheet(worksheet_filename)
        except:
            print('Missing Worksheet Data')

    def open_json(self, filename):
        with open(filename, 'r') as f:
            return json.load(f)

    def open_worksheet(self, filename):
        return pd.read_excel(filename, sheet_name='Samples')

    def get_report_data(self):
        results = {
            'Sequence': self.get_seq_name(),
            'Total Reads': self.get_total_reads(),
            'Lineage': self.get_pango_lineage(),
            'Clade': self.get_clade(),
            'Strain': self.get_strain(),
            'Trim Percent': self.get_trimmed_reads(),
            'Copies per mL': self.get_concentrations(),
            'Base Quality': '99',#self.get_base_quality(),
            'Contamination Rate': self.get_contamination(),
            'Map Rate': self.get_mapped_reads(),
            'Duplicate Read Rate': self.get_duplicate_read_rate(),
            'Median Coverage': self.get_median_coverage(),
            'Average Coverage': self.get_average_coverage(),
            'Positions with Low Coverage Percent': self.get_bases_with_low_coverage(),
            'Number of Substitutions': self.get_substitutions_number(),
            'Number of Deletions': self.get_deletions_number(),
            'Number of Insertions': self.get_insertions_number(),
            'Number of Missing Bases': self.get_missing_bases(),
            'Substitutions': self.get_substitutions(),
            'Deletions': self.get_deletions(),
            'Insertions': self.get_insertions()
        }
        return results

    def validate_worksheet(self):
        pass

    def get_seq_name(self):
        return self.summary['summary_statistics']['table_1'][0]['sample_name'].split('_')[0]

    def get_total_reads(self):
        '''Returns the total reads'''
        return self.summary['summary_statistics']['table_1'][0]['reads_number_of']

    def get_trimmed_reads(self):
        '''Returns the % of reads left after trimming'''
        return self.trim['trim_summary']['table_1'][0]['reads_after_trim_percent']/100

    # def get_base_quality(self):
    #     '''Returns % of base reads with a Q score > 30'''
    #     return self.qc['quality_distribution']['quality_distribution_per_sequence']['table_1'][1]['id_30_40']

    def get_mapped_reads(self):
        '''Returns % of reads post-trimm that map to the SARS-CoV-2 (wuhan-1) genome'''
        #todo get # of reads that map to Wuhan-1 and divide by get_trimmed_reads output
        return self.map_summary['reads_summary']['table_1'][1]['mapped_reads_percent']/100

    def get_contamination(self):
        '''Returns % of reads that map to HG38 which will be the contamination %'''
        #todo get # of reads that map to HG38 and divide by get_trimmed_reads output
        return self.map_summary['reads_summary']['table_1'][0]['mapped_reads_percent']/100

    def get_duplicate_read_rate(self):
        '''Gets the # of duplicate reads'''
        return self.duplicates['table_1'][0]['duplicates_percent']

    def get_median_coverage(self):
        '''Returns the median coverage level'''
        return self.qc_targeted['summary']['table_1'][0]['median_coverage']

    def get_average_coverage(self):
        '''Returns the average coverage level'''
        return self.qc_targeted['summary']['table_1'][0]['avg._coverage']

    def get_bases_with_low_coverage(self):
        '''Returns the % of base positions with low coverage, i.e. <Q30'''
        return self.qc_targeted['summary']['table_1'][0]['length_of_target_region_positions_with_low_coverage_percent']

    def get_substitutions(self):
        '''Returns the base substitutions'''
        subs = self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['substitutions'].values
        if type(subs) == str:
            return subs.split(',')
        else:
            return subs

    def get_deletions(self):
        '''Returns the base deletions'''
        dels =self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['deletions'].values
        if type(dels) == str:
            return dels.split(',')
        else:
            return dels

    def get_insertions(self):
        '''Returns the base insertions'''
        ins = self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['insertions'].values
        if type(ins) == str:
            return ins.split(',')
        else:
            return ins

    def get_substitutions_number(self):
        '''Returns the # of substitutions'''
        return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['totalSubstitutions'].values[0]

    def get_deletions_number(self):
        '''Returns the # of substitutions'''
        return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['totalDeletions'].values[0]

    def get_insertions_number(self):
        '''Returns the # of insertions'''
        return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['totalInsertions'].values[0]

    def get_missing_bases(self):
        '''Returns the missing bases in the consensus sequence'''
        return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['totalMissing'].values[0]

    def get_concentrations(self):
        '''Returns the calculated concentration if available'''
        try:
            return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['Calculated Concentration'].values[0]
        except:
            return 'No Data'

    def get_pango_lineage(self):
        '''Returns the pango lineage'''
        try:
            return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['Nextclade_pango'].values[0]
        except:
            return 'No Data'

    def get_clade(self):
        '''Returns Nextclade clade designation'''
        try:
            return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['clade'].values[0]
        except:
            return 'No Data'

    def get_strain(self):
        '''Returns Known strain info'''
        try:
            return self.sample_summary.loc[self.sample_summary['Sequencing ID'] == self.get_seq_name()]['Strain'].values[0]
        except:
            return 'No Data'

worksheet_location =  r"..\Illumina-COVIDSeq\Data\NGS220809.1KM\NGS220809.1KM_worksheet.xlsx"
folder_location = r"..\Illumina-COVIDSeq\Data\NGS220809.1KM\Reports"

def collect_all_data(folder):
    df = pd.DataFrame()
    for dir, dirnames, files in os.walk(folder):
        for file in files:
            p = os.path.join(folder, file)
            data = pd.DataFrame([covidSeqReport(p, worksheet_location).get_report_data()])
            df = df.append(data, ignore_index=True)
    return df

df = collect_all_data(folder_location)
df = df.assign(Correct = df.Strain == df.Lineage)
df.to_csv(r'..\Illumina-COVIDSeq\Data\NGS220809.1KM\script_output\all_data.csv', index=False)

sensitivity = df.dropna(subset=['Strain']).groupby(['Copies per mL']).sum().reset_index()[['Copies per mL', 'Correct']]
sensitivity = sensitivity.assign(Total = 9)
sensitivity = sensitivity.assign(CI_Lower = lambda x: proportion_confint(x.Correct, x.Total, 0.05, 'wilson')[0])
sensitivity = sensitivity.assign(CI_Upper = lambda x: proportion_confint(x.Correct, x.Total, 0.05, 'wilson')[1])
sensitivity = sensitivity.assign(Percent_Detected = sensitivity.Correct/9)
sensitivity.to_csv(r"C:\DevLab\RND\NGS\Validation\Illumina-COVIDSeq\Data\NGS220809.1KM\script_output\sensitivity.csv", index=False)

cc = df.iloc[:20][['Sequence', 'Lineage', 'Strain', 'Median Coverage', 'Number of Missing Bases', 'Contamination Rate']].rename(columns={'Strain':'Expected'})
cc.to_csv(r'..\Illumina-COVIDSeq\Data\NGS220809.1KM\script_output\clinical_correlation.csv', index=False)

print("PPA CI" + str(proportion_confint(9, 9, 0.05, 'wilson')))
print("NPA CI" + str(proportion_confint(11, 11, 0.05, 'wilson')))
