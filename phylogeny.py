#! /usr/bin/python
import re,os,random,os.path,time
import pandas as pd
#from gff2longest_CDS import GFF2LONGEST_CDS as seq
from fasta import fa_dic
from orthologs_groups import sonic_group_df,ORTHOMCL,busco_group_df
from multiprocessing import Pool

class PHYLOGENY(object):
    def __init__(self,orthocluster_out,species_number,max_gene_number,software_name,raw_database,output_dir,TSCG,PSCG,LCG,OutGroup,DataType,min_PSCG_SPnumber,trim):
        self.orthocluster_out = orthocluster_out #同源基因聚类输出
        self.species_number = species_number # 总共的物种数目
        self.max_gene_number = max_gene_number # 每个family的最大基因个数
        self.software_name = software_name  # 同源聚类的方法
        self.raw_database = raw_database # 原始clean序列的目录，用于提取基因序列
        self.output_dir = output_dir # 输出目录
        self.TSCG = TSCG # 用于提取全部单拷贝基因的物种列表，eg['Atr','Ath']
        self.PSCG = PSCG # 用于提取部分单拷贝基因的物种列表
        self.LCG = LCG # 用于提取低拷贝基因的物种列表
        self.InterGroup = TSCG + PSCG + LCG
        self.OutGroup = OutGroup # 外类群列表
        self.DataType = DataType # 建树数据类型
        self.trim = trim # 用于pxclsq 中修建的比例和trimal
        self.min_PSCG_SPnumber = min_PSCG_SPnumber # 最少内类群物种数目(最大为总共的物种数目-1)
        if self.DataType == 'CDS' or self.DataType == 'Codon12': # iqtree 中数据类型参数
            self.st = 'DNA'
        elif self.DataType == 'pep':
            self.st = 'AA'
        if self.max_gene_number == 1: # 单拷贝还是多拷贝
            self.Type = 'single'
        else:
            self.Type = 'low'
	
    def core_orthogroups(self): # 从orthologs_groups模块中导入标准化后的基因家族聚类数据框，返回过滤后的数据框
        if self.software_name == 'orthomcl':
            all_group_df = ORTHOMCL().all_group_df(self.orthocluster_out)
        elif self.software_name == 'sonicparanoid': # return a DataFrame of single copy gene.
            all_group_df = sonic_group_df(self.orthocluster_out,self.species_number)
        elif self.software_name == 'busco':
            all_group_df = busco_group_df(self.orthocluster_out)
        all_group_df = all_group_df[self.OutGroup + self.InterGroup]
        d = {}
        for groupID,row in all_group_df.iterrows():
            groupID = 'group_' + str(groupID)
            outgroup = []
            for q in self.OutGroup:
                if len(row[q]) != 1: #将超过1个同源基因的外类群定位'*'
                    row[q] = '*'
                elif row[q] != '*': #即长度不等于1又不为'*'的，表示只含有1个同源基因，遂添加到列表，不符合的外类群列表为空
                    outgroup.append(q)
            if len(outgroup) > 1: #如果外类群个数超过1，则只选用第一个物种的基因作为外类群，目前只支持两个外类群物种
                row[outgroup[1]] = '*'
            TSCG = []
            PSCG = []
            LCG = []
            for i in self.TSCG:
                if row[i] != '*' and len(row[i]) == 1:
                    TSCG.append(i)
                else:
                    row[i] = '*'
            for i in self.PSCG:
                if row[i] != '*' and len(row[i]) == 1:
                    PSCG.append(i)
                else:
                    row[i] = '*'
            for i in self.LCG:
                if row[i] != '*' and len(row[i]) <= self.max_gene_number:
                    LCG.append(i)
                else:
                    row[i] = '*'
            if len(TSCG) == len(self.TSCG) and outgroup: #严格单拷贝：所有TSCG物种都是单拷贝
                if len(PSCG) >= self.min_PSCG_SPnumber: #部分单拷贝：该部分物种有一定比例的单拷贝即可，self.min_PSCG_SPnumber等于内类群的总个数，则所有物种都是单拷贝(严格单拷贝基因集)
                    if len(LCG) == len(self.LCG):
                        d[groupID] = row
        core = pd.DataFrame(d).T #.T为横纵坐标转置
        for abb,row in core.iteritems():
            print(abb+'\t'+str(len([m for m in row if m != '*'])/len(row)),len(row))
        print('output the groups csv')
        core.to_csv(self.software_name+'_'+self.Type+'.csv',sep='\t')
        return core # 返回筛选后的数据框
    
    def extract_tree(self):
        #core_orthogroups = self.core_orthogroups()
        df = pd.read_csv('/data/01/user119/project/jsl/04_phylogeny/orthomcl/orthomcl_single.csv',sep='\t',index_col=0)
        for groupID in df.index.values:
            tree_path = os.path.join(self.output_dir,'align',groupID,groupID+'_'+self.DataType+'.aln-trimal.treefile')
            if os.path.exists(tree_path):
                print(tree_path)

    def extraction_singlecopy_cds(self,input_cds_fa,abb):
        fa_seq = fa_dic(input_cds_fa)
        singlecopy_geneID = self.singlecopy_geneID(orthomcl_singlecopy)
        pass       

    def build_database_1(self): # 生成每个基因的CDS与pep文件
        core_orthogroups = self.core_orthogroups()
        if not os.path.exists(self.output_dir+'/align'):
            os.mkdir(self.output_dir+'/align')
        for groupID,row in core_orthogroups.iterrows():
            print(groupID)
            path = os.path.join(self.output_dir,'align',groupID)
            CDS_path = os.path.join(self.output_dir,'align',groupID,groupID+'_CDS')
            pep_path = os.path.join(self.output_dir,'align',groupID,groupID+'_pep')
            if not os.path.exists(CDS_path) or not os.path.exists(pep_path): #如果没有生产pep文件，则新建
                if not os.path.exists(path):
                    os.mkdir(path)
            else:
                continue # 如果有，则跳过
            f_cds = open(CDS_path,'w')
            f_pep = open(pep_path,'w')
            for geneID_list in row:
                if geneID_list == '*':continue
                for index,geneID in enumerate(geneID_list):
                    spID = geneID[:3]
                    if spID in self.OutGroup:
                        spID = 'Out'
                    cds_seq = fa_dic(os.path.join(self.raw_database,geneID[:3]+'_cds.fa'))[geneID]
                    pep_seq = fa_dic(os.path.join(self.raw_database,geneID[:3]+'_pep.fa'))[geneID]
                    f_cds.write('>'+spID+'_'+str(index+1)+'\n'+cds_seq+'\n')
                    f_pep.write('>'+spID+'_'+str(index+1)+'\n'+pep_seq+'\n')
            f_cds.close()
            f_pep.close()
        
    def build_database_multiprocess_1(self,row): # 生成每个基因的CDS与pep文件
        # row 为 (index,series)
        if not os.path.exists(self.output_dir+'/align'):
            os.mkdir(self.output_dir+'/align')
        groupID = row[0]
        path = os.path.join(self.output_dir,'align',groupID)
        CDS_path = os.path.join(self.output_dir,'align',groupID,groupID+'_CDS')
        pep_path = os.path.join(self.output_dir,'align',groupID,groupID+'_pep')
        if not os.path.exists(CDS_path) or not os.path.exists(pep_path): #如果没有生成pep文件，则新建
            if not os.path.exists(path):
                os.mkdir(path)
        else:
            exit(0) # 如果有，则跳过
        f_cds = open(CDS_path,'w')
        f_pep = open(pep_path,'w')
        for geneID_list in row[1]:
            if geneID_list == '*':continue
            for index,geneID in enumerate(geneID_list):
                spID = geneID[:3]
                if spID in self.OutGroup:
                    spID = 'Out'
                cds_seq = fa_dic(os.path.join(self.raw_database,geneID[:3]+'_cds.fa'))[geneID]
                pep_seq = fa_dic(os.path.join(self.raw_database,geneID[:3]+'_pep.fa'))[geneID]
                f_cds.write('>'+spID+'_'+str(index+1)+'\n'+cds_seq+'\n')
                f_pep.write('>'+spID+'_'+str(index+1)+'\n'+pep_seq+'\n')
        f_cds.close()
        f_pep.close()
        print(groupID)

    def singleCOPY_align_2(self): # 单个基因序列比对
        with open(os.path.join(self.output_dir,'align_%s_%s.sh' % (self.DataType,self.trim)),'w') as f:
            for groupID in os.listdir(self.output_dir+'/align'):
                path = os.path.join(self.output_dir+'/align',groupID)
                if not os.path.exists(os.path.join(path,groupID+'_CDS.aln')) or not os.path.exists(os.path.join(path,groupID+'_pep.aln')):
                    #f.write('cd %s && mafft --auto --quiet %s_pep > %s_pep.aln && pal2nal.pl %s_pep.aln %s_CDS -output fasta > %s_CDS.aln && trimal -in %s_pep.aln  -out %s_pep.aln-trimal -automated1 && trimal -in %s_CDS.aln -out %s_CDS.aln-trimal -automated1\n'  % (path,groupID,groupID,groupID,groupID,groupID,groupID,groupID,groupID,groupID))
                    f.write('cd %s && mafft --auto --quiet %s_pep > %s_pep.aln && pal2nal.pl %s_pep.aln %s_CDS -output fasta > %s_CDS.aln && pxclsq -a -p 0.2 -s %s_pep.aln  -o %s_pep.aln-0.2cln && pxclsq -a -p 0.2 -s %s_CDS.aln -o %s_CDS.aln-0.2cln\n' % (path,groupID,groupID,groupID,groupID,groupID,groupID,groupID,groupID,groupID))

    def aligned_filter_3(self): # 过滤掉非保守位点
        for groupID in os.listdir(self.output_dir+'/align'):
            print(groupID)
            aligned_fa = os.path.join(self.output_dir+'/align',groupID,groupID+'_%s_aligned' % self.DataType)
            d = []
            fa_seq = fa_dic(aligned_fa)
            t = {m:list(n) for m,n in fa_seq.items()}
            matrix = pd.DataFrame(t)
            for index,row in matrix.iterrows():
                if row.tolist().count('-') <= 3:
                    d.append(row)
            if len(d) < 300:continue
            f = open(aligned_fa + '_filtered','w')
            matrix_filtered = pd.DataFrame(d)
            for index,row in matrix_filtered.iteritems():
                f.write('>'+index+'\n'+''.join(row.tolist())+'\n')
            f.close()

    def CDS2Codon12(self): # 根据CDS提取codon12
        for groupID in os.listdir(self.output_dir+'/align'):
            path = os.path.join(self.output_dir,'align',groupID)
            if not os.path.exists(os.path.join(path,groupID+'_CDS.aln-0.2cln')):
                continue
            f_codon12 = open(os.path.join(path,groupID+'_Codon12.aln-0.2cln'),'w')
            Fa_dic = fa_dic(os.path.join(path,groupID+'_CDS.aln-0.2cln'))
            for geneID,seq in Fa_dic.items():
                d = []
                for n in range(0,len(seq),3):
                    d.append(seq[n:n+2])
                f_codon12.write('>'+geneID+'\n'+''.join(d)+'\n')
            f_codon12.close()
            print(groupID)

    def genetree_iqtree_4(self): # 生成构建基因树sh
        f = open(os.path.join(self.output_dir,'tree_%s_%s.sh' % (self.DataType,self.trim)),'w')
        for groupID in os.listdir(self.output_dir+'/align'):
            path = os.path.join(self.output_dir,'align',groupID)
            if os.path.exists('%s/%s.aln-%s.treefile' % (path,groupID+'_'+self.DataType,'0.2cln')):continue
            if os.path.exists('%s/%s.aln-%s' % (path,groupID+'_'+self.DataType,"0.2cln")):
                f.write('cd %s && iqtree -s %s.aln-%s  -pre  %s.aln-%s -st %s  -nt  AUTO  -bb 1000 -redo -m MFP -quiet\n' % (path,groupID+'_'+self.DataType,self.trim,groupID+'_'+self.DataType,self.trim,self.st))
        f.close()

    def generate_concatenation_5(self): # 将序列连接成超级矩阵
        core_orthogroups = self.core_orthogroups()
        print(len(core_orthogroups))
        d = {i[:3]:[] for i in self.InterGroup+['Out']}
        for groupID,row in core_orthogroups.iterrows():
            path = os.path.join(self.output_dir,'align',groupID,groupID+'_%s.aln-0.2cln' % (self.DataType)) 
            if os.path.exists(path):
                t = fa_dic(path)
                spID_current = set([i[:3] for i in t.keys()])
                for x in set(d)-spID_current:
                    d[x].append('-'*len(t[random.choice(list(t.keys()))]))
                for m,n in t.items():
                    d[m[:3]].append(n)
        with open(self.output_dir+'/%s_concatenation.fa.aln-0.2cln' % (self.DataType),'w') as f:
            for i in d:
                f.write('>'+i+'\n'+''.join(d[i])+'\n')

def generate_concatenation(suffix,aln_dir,out_name):
    d = {}
    for i in os.listdir(aln_dir):
        if i.endswith('CDS.aln-0.2cln'):
            fa = fa_dic(os.path.join(aln_dir,i))
            d[i] = fa
    df = pd.DataFrame(d).fillna('*')
    current_len = 0
    with open(out_name,'w') as f:
        for abb,row in df.iterrows():
            for groupID,s in row.items():
                gene_length = len(list(fa_dic(os.path.join(aln_dir,groupID)).values())[0])
                if s == '*':
                    gap = '-'*gene_length
                    row[groupID] = gap
                print('DNA, '+groupID+'='+str(current_len+1)+'-'+str(current_len+gene_length))
                current_len += gene_length
            seq = ''.join(row.tolist())
            f.write('>'+abb+'\n'+seq+'\n')

def write_seq_from_tree(self,output_dir,tree_dir,pep_file,CDS_file): # 提取树文件中基因的序列
    pep_dic = fa_dic(pep_file)
    CDS_dic = fa_dic(CDS_file)
    if not os.path.exists(output_dir+'/align'):
        os.mkdir(output_dir+'/align')
    for i in os.listdir(tree_dir):
        cluster_name = i.split('.')[0]
        path = os.path.join(output_dir,'align',cluster_name)
        if not os.path.exists(path):
            os.mkdir(path)
        print(cluster_name)
        #f_list = open(cluster_name+'_gene_list','w')
        treestr = open(os.path.join(tree_dir,i)).read()
        all_gene = re.findall('(\w{3}@.*?):',treestr)
        with open(os.path.join(path,cluster_name+'_pep'),'w') as f_pep, open(os.path.join(path,cluster_name+'_CDS'),'w') as f_CDS:
            for geneID in all_gene:
                f_pep.write('>'+geneID[:3]+'\n'+pep_dic[geneID]+'\n')
                f_CDS.write('>'+geneID[:3]+'\n'+CDS_dic[geneID]+'\n')
            #f_list.write(geneID+'\n')
        #f_list.close()

def extract_seq_for_treeshrink(clusterID,outdir): #提取 treeshink 删除基因后的对齐序列
    for clusterID,row in pd.read_csv(clusterID,sep='\t',header=None).set_index(0).fillna('*').iterrows():
        print(clusterID)
        cds_path = os.path.join('/data/01/user119/project/jsl/04_phylogeny/orthomcl/Gbi_Pab_single/align/%s/%s_CDS.aln-0.2cln' % (clusterID,clusterID))
        fa = fa_dic(cds_path)
        for i in row:
            if i != '*':
                del fa[i]
        with open(os.path.join(outdir,clusterID + '_CDS_shrunk_0.2cln.fasta'), 'w') as f:
            for geneID,seq in fa.items():
                f.write('>'+geneID+'\n'+seq+'\n')


if __name__ == '__main__':
    #raw_database = ''
    #output_dir = ''
    #orthomcl_out = ''
    #output_dir = ''
    #a = ORTHOMCL()
	#d = PHYLOGENY(orthomcl_out,15,1,'orthomcl',raw_database,output_dir,'Efe',['Gbi','Pab'],'CDS',12)
	#a.core_orthogroups(orthomcl_out,14,1)
	#d.build_database_1()
	#d.singleCOPY_align_2()
	#d.aligned_filter_3()
	#d.genetree_iqtree_4()
	#d.generate_concatenation_5()
	#d.CDS2Codon12()
	#d.extraction_singlecopy_cds()
    #generate_concatenation()