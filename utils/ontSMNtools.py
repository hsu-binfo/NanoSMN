import re
import os
import io
import glob
import gzip
import errno
import jinja2
import requests, sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce

class ontSMNtools:
  def __init__(self, path):
    self.__mapping = None
    self.__path = path
    self.get_SMN_mapping_table()

    if not os.path.isdir(path):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)

  def split_reads_by_type(self, sample):
    all_gaf = glob.glob(os.path.join(self.__path, sample, "*.gaf"))
    bin = {}
    for gaf in all_gaf:
      read_types = self.__get_read_types_from_gaf(gaf)
      for read_item in read_types:
        if read_item['type'] not in bin:
          bin[read_item['type']] = []
        bin[read_item['type']].append(read_item['id'])
    self.__filter_fastq(sample, bin)

  def plot_cnv(self, sample):
    all_gaf = glob.glob(os.path.join(self.__path, sample, "*.gaf"))
    for gaf in all_gaf:
      # sample = self.__get_sample_name(gaf)
      read_types = self.__get_read_types_from_gaf(gaf)
      df_stat = pd.DataFrame(read_types).groupby(['type', 'node_index']).agg(
        SMN1 = ('SMN1', sum),
        SMN2 = ('SMN2', sum)
      ).reset_index()
      
      df_stat['total'] = df_stat['SMN1'] + df_stat['SMN2']
      df_stat['SMN1_percentage'] = df_stat['SMN1'] / df_stat['total']
      df_stat['SMN2_percentage'] = df_stat['SMN2'] / df_stat['total']
      df_stat['ratio'] = df_stat['SMN2'] / df_stat['SMN1']
      df_stat['ratio_log2'] = np.log2(df_stat['ratio'])
      df_stat['correct_ratio_log2'] = df_stat['ratio_log2'].apply(self.__correct_ratio)

      (SMN1, SMN2, SMN_ratio) = self.__get_SMN_ratio(read_types)
      self.__generate_plot(sample, df_stat, SMN_ratio)

  def get_SMN_mapping_table(self):
    SMN2="1+,5693+,5694+,6104+,6105+,7477+,8512+,8513+,9056+,9057+,10221+,10222+,11956+,11957+,12097+,12098+,12101+,12102+,12239+,12240+,12251+,12252+,12405+,12406+,12879+,12880+,12951+,12956+,13011+,13012+,13312+,13313+,13384+,13385+,13819+,13820+,13869+,13870+,14018+,14019+,14134+,14135+,14600+,14601+"
    SMN1="1+,5694+,6103+,6105+,7475+,7477+,8511+,8513+,9055+,9057+,10220+,10222+,11955+,11957+,12096+,12098+,12100+,12102+,12238+,12240+,12250+,12252+,12404+,12406+,12878+,12880+,12956+,13010+,13012+,13311+,13313+,13383+,13385+,13818+,13820+,13868+,13870+,14017+,14019+,14133+,14135+,14599+,14601+,14872+"
    SMN1_nodes = self.__split_path(SMN1)
    SMN2_nodes = self.__split_path(SMN2)

    # select the nodes we are interested in
    df_SMN1 = pd.DataFrame(SMN1_nodes, columns=['SMN1'])
    df_SMN1 = df_SMN1[(df_SMN1['SMN1'] >= 10220) & (df_SMN1['SMN1'] <= 14600)]
    df_SMN1.reset_index(inplace = True, drop = True)
    self.__df_SMN1 = df_SMN1
    df_SMN2 = pd.DataFrame(SMN2_nodes, columns=['SMN2'])
    df_SMN2 = df_SMN2[(df_SMN2['SMN2'] >= 10220) & (df_SMN2['SMN2'] <= 14600)]
    df_SMN2.reset_index(inplace = True, drop = True)
    self.__df_SMN2 = df_SMN2

    df_merge = df_SMN1.merge(df_SMN2, how="outer", left_on="SMN1", right_on="SMN2")
    df_merge['max'] = df_merge[["SMN1", "SMN2"]].max(axis=1)
    df_merge = df_merge.sort_values('max')

    mapping = self.__get_node_mapping_table(df_merge)
    self.__mapping = mapping
    return self.__mapping
  
  def variant_annotation(self, sample, type):
    vcf_file = os.path.join(self.__path, sample, type, "variant_calls.final.vcf.gz")
    with gzip.open(vcf_file, 'r') as f:
      lines = [l.decode('utf-8') for l in f if not l.decode('utf-8').startswith('##')]

    vcf = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

    if vcf.shape[0]:
      vcf['CHROM_POS'] = vcf.apply(self.__map_ref_position, axis=1)
      self.__variant_annotation(sample, type, vcf)
    else:
    # no variant called --> generate final file for validation
      sample_path = os.path.join(self.__path, sample)
      output_file = os.path.join(sample_path, "{}.variant_calls.final.vep".format(type))
      f = open(output_file, "w")
      f.close()

  def generate_report(self, sample):
    template_path = os.path.dirname(__file__)
    template_loader = jinja2.FileSystemLoader(searchpath=template_path, encoding='utf8')
    template_env = jinja2.Environment(loader=template_loader)
    TEMPLATE_FILE = "report.html"
    template = template_env.get_template(TEMPLATE_FILE)

    stat_file = os.path.join(self.__path, sample, sample+".stat")
    (raw_len, raw_qual) = self.__parse_stat(stat_file)

    stat_file = os.path.join(self.__path, sample, sample+".filtered.stat")
    (filter_len, filter_qual) = self.__parse_stat(stat_file)
    cnv_plot = sample+".cnv.png"
    
    gaf = os.path.join(self.__path, sample, sample+".aln.gaf")
    read_types = self.__get_read_types_from_gaf(gaf)
    (SMN1, SMN2, SMN_ratio) = self.__get_SMN_ratio(read_types)

    SMN1_vep = os.path.join(self.__path, sample, "SMN1.variant_calls.final.vep")
    SMN1_vep_data = pd.read_csv(SMN1_vep, sep="\t").to_html(index=False, escape=False, border=0, justify='left', classes='table table-striped table-hover')
    SMN2_vep = os.path.join(self.__path, sample, "SMN2.variant_calls.final.vep")
    SMN2_vep_data = pd.read_csv(SMN2_vep, sep="\t").to_html(index=False, escape=False, border=0, justify='left', classes='table table-striped table-hover')

    content = template.render(
      sample=sample,
      raw_len=raw_len,
      raw_qual=raw_qual,
      filter_len=filter_len,
      filter_qual=filter_qual,
      SMN1=SMN1,
      SMN2=SMN2,
      SMN_ratio=round(SMN_ratio, 2),
      cnv_plot=cnv_plot,
      SMN1_vep=SMN1_vep_data,
      SMN2_vep=SMN2_vep_data
    )
    report_file = os.path.join(self.__path, sample, "report.html")
    with open(report_file, "w") as f:
      f.write(content)

  def __split_path(self, path):
    all_nodes = path.split(",")
    node_list = []
    for node in all_nodes:
      m = re.search(r"(\d+)([+-])", node)
      pos = m.group(1)
      strand = m.group(2)
      node_list.append(int(pos))
    return node_list

  def __merge_nodes(self, read_nodes):
    all_nodes = []
    for nodes in read_nodes:
      for node in nodes:
        all_nodes.append(node)
    return list(dict.fromkeys(all_nodes))

  def __get_node_mapping_table(self, df_merge):
    last_SMN1 = None
    last_SMN2 = None
    flag = 0
    mapping = []
    for index, row in df_merge.iterrows():
      if pd.notnull(row['SMN1']) and pd.notnull(row['SMN2']):
        if not flag and (last_SMN1 or last_SMN2):
          mapping.append({
            "SMN1": last_SMN1,
            "SMN2": last_SMN2
          })
        # mapping.append({
        #   "SMN1": row['SMN1'],
        #   "SMN2": row['SMN2']
        # })
        flag = 0
      elif pd.isnull(row['SMN1']) and pd.isnull(last_SMN2):
        mapping.append({
          "SMN1": last_SMN1,
          "SMN2": row['SMN2']
        })
        flag = 1

      last_SMN1 = row['SMN1']
      last_SMN2 = row['SMN2']
    if flag == 0:
      mapping.append({
        "SMN1": last_SMN1,
        "SMN2": last_SMN2
      })
    return mapping

  def __get_read_types_from_gaf(self, gaf):
    read_nodes = {}
    read_types = []
    with open(gaf) as f:
      for line in f:
        columns = line.rstrip().split("\t")
        path = columns[5]
        id = columns[0].split(" ")[0]
        if id not in read_nodes:
          read_nodes[id] = []

        if re.match(">", path):
          nodes = path.split(">")
          nodes.remove('')
        elif re.match("<", path):
          nodes = path.split("<")
          nodes.remove('')
          nodes.reverse()
        read_nodes[id].append(nodes)

    for id in read_nodes:
      all_nodes = self.__merge_nodes(read_nodes[id])
      (score_SMN1, mapping_SMN1) = self.__infer_mapping_score(all_nodes, list(self.__df_SMN1['SMN1']))
      (score_SMN2, mapping_SMN2) = self.__infer_mapping_score(all_nodes, list(self.__df_SMN2['SMN2']))
      if score_SMN1 > score_SMN2:
        type = "SMN1"
      else:
        type = "SMN2"
      node_types = self.__find_node_types(all_nodes)
      for index, node_type in node_types.items():
        if node_type == "SMN1":
          read_types.append({
            "id": id,
            "type": type,
            "node_index": index,
            "SMN1": 1,
            "SMN2": 0
          })
        else:
          read_types.append({
            "id": id,
            "type": type,
            "node_index": index,
            "SMN1": 0,
            "SMN2": 1
          })
    return read_types

  def __get_SMN_ratio(self, read_types):
    SMN1 = pd.DataFrame(read_types).groupby(['type']).nunique()['id']['SMN1']
    SMN2 = pd.DataFrame(read_types).groupby(['type']).nunique()['id']['SMN2']
    SMN_ratio = SMN1/SMN2

    return SMN1, SMN2, SMN_ratio

  def __infer_mapping_score(self, nodes, ref_nodes):
    all_matches = []
    match = []
    score = 0
    last_index = None
    for node in nodes:
      node = int(node)
      if node in ref_nodes:
        index = ref_nodes.index(node)
        if last_index is None:
          match = []
          # match.append(node)
        elif index == last_index + 1:
          # match.append(node)
          pass
        else:
          all_matches.append(match)
          match = []
        match.append(node)

        score +=1
        last_index = index
      else:
        continue
    if len(match) > 0:
      all_matches.append(match)
    return score, all_matches

  def __find_node_types(self, all_nodes):
    result = {}
    mapping = self.__mapping
    for m in mapping:
      index = mapping.index(m)
      matched = 0
      for type in m.keys():
        if not np.isnan(m[type]) and str(int(m[type])) in all_nodes:
          result[index] = type
          matched = 1
      if not matched:
        for type, pos in m.items():
          if np.isnan(pos):
            result[index] = type
    return result

  def __correct_ratio(self, x, limit=1.5):
    if x > limit:
      return limit
    elif x < -limit:
      return -limit
    return x

  def __generate_plot(self, sample, df_draw, SMN_ratio):
    sample_path = os.path.join(self.__path, sample)
    if not os.path.isdir(sample_path):
      os.mkdir(sample_path)

    sns.set_style("ticks")
    plt.rcParams["figure.figsize"] = (12,4)
    plt.rcParams.update({'font.size': 14})

    mapping = self.__mapping

    # 計算SMN1/SMN2比例
    df_ratio = df_draw.groupby('node_index').agg({'total': lambda T: reduce(np.divide, T)})
    df_ratio.rename(columns = {'total':'SMN1/SMN2'}, inplace = True)
    df_ratio['SMN1/SMN2_log'] = np.log2(df_ratio['SMN1/SMN2'])
    df_ratio['SMN1/SMN2_log'] = df_ratio['SMN1/SMN2_log'].apply(lambda x: self.__correct_ratio(x, 2.5))
    SMN_ratio_log = self.__correct_ratio(np.log2(SMN_ratio), 2.5)
    display_ratio_string = 'SMN1/SMN2='+str(np.round(SMN_ratio, 2))
    if SMN_ratio_log == -2.5:
      display_ratio_string = 'SMN1 Deletion'
      df_draw.loc[df_draw['type'] == 'SMN1', 'correct_ratio_log2'] = -1.5
    elif SMN_ratio_log == 2.5:
      display_ratio_string = 'SMN2 Deletion'
      df_draw.loc[df_draw['type'] == 'SMN2', 'correct_ratio_log2'] = 1.5

    g = sns.lineplot(data=df_draw, x="node_index", y="correct_ratio_log2", hue="type", marker="o", markersize=10)
    g.set_xticks(range(len(mapping)))
    g.set_xticklabels([
      # '5693:-/T',
      # '6103/6104:G/A',
      # '7475:TG/-',
      # '8511/8512:T/C',
      # '9055/9056:A/G',
      '10220/10221:A/G',
      '11955/11956:T/C',
      '12096/12097:G/A',
      '12100/12101:T/C',
      '12238/12239:G/A',
      '12250/12251:T/C',
      '12404/12405:G/A',
      '12878/12879:G/A',
      '12951:----/GCAAG',
      '13010/13011:A/C',
      '13311/13312:G/A',
      '13383/13384:T/C',
      '13818/13819:G/A',
      '13868/13869:C/T',
      '14017/14018:A/G',
      '14133/14134:A/G',
      '14599/14600:G/A',
      # '14872:G/-'
      ], rotation = 90)

    g.set_title(sample, y=1.1, weight='bold')
    g.set_ylim([-2.8, 2.8])
    g.yaxis.set_ticks(np.arange(-2, 3, 1))
    g.legend(title=None)
    g.xaxis.set_label_text(" ")
    g.yaxis.set_label_text('log2(Ratio)')
    plt.axhline(y=0, color='silver', linestyle='--', linewidth=1.3)
    plt.axhline(y=1, color='silver', linestyle='--', linewidth=1.3)
    plt.axhline(y=-1, color='silver', linestyle='--', linewidth=1.3)

    df_ratio.drop(columns=['SMN1/SMN2'], inplace=True)
    df_ratio.rename(columns={'SMN1/SMN2_log': 'SMN1/SMN2'}, inplace=True)
    sns.lineplot(data=df_ratio, palette=['g'])
    plt.text(len(mapping)-4, SMN_ratio_log+0.2, display_ratio_string, fontsize=14)

    # 設定ledgend位置
    g.legend(bbox_to_anchor=(1.01, 1.36))

    output_file = os.path.join(sample_path, "{}.cnv.png".format(sample))
    plt.savefig(output_file, bbox_inches='tight')
    plt.clf()

  def __filter_fastq(self, sample, bin):
    sample_path = os.path.join(self.__path, sample)
    if not os.path.isdir(sample_path):
      os.mkdir(sample_path)
    for type in bin.keys():
      bin[type] = list(set(bin[type]))
      output_file = os.path.join(sample_path, "{}_reads.lst".format(type))
      with open(output_file, "w") as f:
        for id in bin[type]:
          f.write(id+"\n")

  def __map_ref_position(self, row):
    start = None
    if row['CHROM'] == 'SMN1':
      start = 70938100
    elif row['CHROM'] == 'SMN2':
      start = 70062676
    return start+row['POS']-1

  def __run_vep(self, ext, return_data, allele_freq):
    server = "https://rest.ensembl.org"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
  
    vep = r.json()[0]
    # SMN1: ENST00000380707
    # SMN2: ENST00000380743
    
    target_transcripts = ['ENST00000380707', 'ENST00000380743']
    target_consequence_obj = next(
      (m for m in vep['transcript_consequences'] if m.get('transcript_id') in target_transcripts),
      None
    )
    return_data += [
      vep['allele_string'],
      allele_freq,
      ",".join(target_consequence_obj['consequence_terms'])
    ]

    cds_change = None
    protein_change = None
    if 'hgvsc' in target_consequence_obj:
      cds_change = target_consequence_obj['hgvsc'].split(":")[1]
    if 'hgvsp' in target_consequence_obj:
      protein_change = target_consequence_obj['hgvsp'].split(":")[1]
    if 'intron' in target_consequence_obj:
      intron_number = target_consequence_obj['intron'].split("/")[0]
      info = 'ISV{} {}'.format(intron_number, cds_change)
    else:
      exon_number = target_consequence_obj['exon'].split("/")[0]
      info = "EXON {} {};{}".format(exon_number, cds_change, protein_change)
    return_data += [info]

    existing_variations = []
    clinical_significants = []
    if 'colocated_variants' in vep:
      for v in vep['colocated_variants']:
        ## Existing_variation
        if re.match('^rs', v['id']):
          existing_variations.append('<a href="https://www.ncbi.nlm.nih.gov/snp/{}" target="_blank">{}</a>'.format(v['id'], v['id']))
        elif re.match('^COSV', v['id']):
          existing_variations.append('<a href="https://cancer.sanger.ac.uk/cosmic/search?q={}" target="_blank">{}</a>'.format(v['id'], v['id']))
        else:
          existing_variations.append(v['id'])

        ## ClinVar
        clinvar_link = []
        if 'var_synonyms' in v:
          clinvar_list = v['var_synonyms'].get('ClinVar')
          for clinvar_id in clinvar_list:
            if re.match('^RCV', clinvar_id):
              clinvar_link.append('<a href="https://www.ncbi.nlm.nih.gov/clinvar/{}/" target="_blank">{}</a>'.format(clinvar_id, clinvar_id))

          clin_sig = "{} ({})".format(
            '/'.join(clinvar_link),
            '/'.join(v.get('clin_sig'))
          )
          clinical_significants.append(clin_sig)
    return_data += [",".join(existing_variations), ",".join(clinical_significants)]

      
    extra = "IMPACT={};STRAND={}".format(target_consequence_obj['impact'], target_consequence_obj['strand'])

    return_data.append(extra)
    return return_data

  def __variant_annotation(self, sample, type, vcf):
    sample_path = os.path.join(self.__path, sample)
    output_file = os.path.join(sample_path, "{}.variant_calls.final.vep".format(type))
    f = open(output_file, "w")

    headers = [
      'Location', 'Allele', 'Genotype(GT:DP:FQ)', 'Consequence', 'Information', 'Existing_variation', 'ClinVar', 'Extra'
    ]
    f.write("\t".join(headers)+"\n")
    for item in vcf.itertuples():
      start = item.CHROM_POS
      end = item.CHROM_POS-1+len(item.REF)
      alts = item.ALT.split(',')
      freq = item.SAMPLE

      for alt in alts:
        ext = "/vep/human/region/5:{}:{}/{}?numbers=1&hgvs=1".format(start, end, alt)
      
        return_data = ["5:{}".format(start)]
        return_data = self.__run_vep(ext, return_data, freq)

        f.write("\t".join(return_data)+"\n")
    f.close()

  def __parse_stat(self, file):
    data = {}
    with open(file, "r") as f:
      lines = f.readlines()
      general_stat = lines[1:9]
      qual_stat = lines[10:15]

      for stat in general_stat:
        stat = stat.rstrip()
        (key, value) = stat.split(":")
        data[key] = value.replace(' ', '')

      for stat in qual_stat:
        stat = stat.rstrip()
        (key, value) = stat.split(r":")
        value = value.replace("\t", '')
        data[key] = value

    len = [data['Mean read length'], data['Median read length'], data['STDEV read length'], data['Read length N50']]
    qual = [data['Mean read quality'], data['Median read quality'], data['>Q5'], data['>Q7'], data['>Q10'], data['>Q12'], data['>Q15']]

    return len, qual