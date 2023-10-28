import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.gridspec import GridSpec
from networkx.algorithms import bipartite
from networkx.drawing.layout import bipartite_layout
from scipy.optimize import curve_fit
import random
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from bs4 import BeautifulSoup
import requests
import urllib.request
import xml.etree.ElementTree as ET
from pyvis.network import Network
import streamlit.components.v1 as components
import zipfile
import os
import scipy.io
import scipy.stats as stats
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as multi
from matplotlib_venn import venn3
from PIL import Image
import statsmodels.api as sm
import gc
#plt.rcParams['font.family']= 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['font.size'] = 12 
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['figure.dpi'] = 300

#Main
st.set_page_config(layout="wide")
st.title("iTraNet: integrated Trans-omics Network visualization and analysis")
st.write("This website is free and open to all users and there is no login requirement.")

#1, upload omics data
#Slidebar
st.sidebar.subheader('1, Upload omics data')
#while not os.path.isfile('demo.zip'):
if(os.path.isfile('demo.zip')):
    os.remove('demo.zip')
with zipfile.ZipFile('demo.zip', 'x') as csv_zip:
    csv_zip.writestr("Transcriptome data (organ or cell).csv",
                    pd.read_csv("./DEMDEG/WT_Liver_Transcriptome.csv").to_csv(index=False))
    csv_zip.writestr("Proteome data (organ or cell).csv",
                    pd.read_csv("./DEMDEG/OBWTdiff_Proteome.csv").to_csv(index=False))
    csv_zip.writestr("Metabolome data (organ or cell).csv", 
                    pd.read_csv("./DEMDEG/WT_Liver_Metabolome.csv").to_csv(index=False))
    csv_zip.writestr("Metabolome data (blood or medium).csv", 
                    pd.read_csv("./DEMDEG/WT_Blood_Metabolome.csv").to_csv(index=False)) 
    csv_zip.writestr("Background gene.csv", 
                    pd.read_csv("./Database/RefGene.csv").dropna().reset_index(drop=True).to_csv(index=False))        
with open("demo.zip", "rb") as file:
    st.sidebar.download_button(label = "Download demo data",data = file,file_name = "demo.zip")

Tran = st.sidebar.file_uploader("Transcriptome (organ or cell)", type="csv")
Pro = st.sidebar.file_uploader("Proteome (organ or cell)", type="csv")
Meta1 =st.sidebar.file_uploader("Metabolome (organ or cell)", type="csv")
Meta_blood=st.sidebar.file_uploader("Metabolome (blood or medium)", type="csv")
df_ref=st.sidebar.file_uploader("Option: background gene (for TF estimation)", type="csv")
if Tran is not None:
    st.subheader('Uploaded transcriptome data (organ or cell)')
    Tran = pd.read_csv(Tran)
    st.write(Tran.set_index('ENSMUSG'))     
if Pro is not None:
    st.subheader('Uploaded proteome data (organ or cell)')
    Pro = pd.read_csv(Pro)
    st.write(Pro.set_index('ENSMUSP'))      
if Meta1 is not None:
    st.subheader('Uploaded metabolome data (organ or cell)')
    Meta1 =pd.read_csv(Meta1)
    st.write(Meta1.set_index('CPD'))        
if Meta_blood is not None:
    st.subheader('Uploaded metabolome data (blood or medium)')
    Meta_blood =pd.read_csv(Meta_blood)
    st.write(Meta_blood.set_index('CPD'))
if df_ref is not None:
    st.subheader('Uploaded background gene (for TF estimation)')
    df_ref =pd.read_csv(df_ref)
    df_ref=df_ref.dropna().reset_index(drop=True)
    st.write(df_ref)
if df_ref is None:
    df_ref=pd.read_csv("./Database/RefGene.csv")
    df_ref=df_ref.dropna().reset_index(drop=True)
   
st.sidebar.subheader('2, Select an analysis method')
options = ['Home','A, gene regulatory network (including TF, miRNA, and mRNA) (transcriptome)',
        "B, mRNA (protein)-mRNA (protein) interaction (transcriptome or proteome)",
        "C, metabolic network (including enzyme, mRNA, and metabolite) (transcriptome and metabolome (organ or cell))", 
        "D, metabolite exchange network (including transporter, mRNA, and metabolite) (transcriptome, metabolome (organ or cell) and metabolome (blood or medium))"]
selected_option = st.sidebar.radio('*"gene regulatory network" takes more than ten minutes',options, key="2")

st.sidebar.subheader('3, Set parameters')
options = ['cyan','blue','green','red','magenta','yellow','black','white']
Inc_color = st.sidebar.selectbox('Increased color:',options)
options = ['red','blue','green','cyan','magenta','yellow','black','white']
Dec_color = st.sidebar.selectbox('Decreased color:',options)
options1 = ['Liver', 'Muscle','All', 'Pluripotent stem cell', 'Unclassified', 'Blood', 'Embryo',
       'Gonad', 'No description', 'Neural','Embryonic fibroblast', 'Digestive tract', 'Kidney',
       'Cardiovascular', 'Pancreas', 'Lung', 'Spleen', 'Bone','Others', 'Uterus', 'Breast', 
       'Epidermis', 'Prostate', 'Adipocyte','Placenta']
Organ= st.sidebar.selectbox('Organ or cell (used for TF estimation):',options1)
Thre=st.sidebar.slider('Threshold (used for TF estimation):', 50, 500, 50)
FDR=st.sidebar.slider('FDR (used for TF estimation):', 0.0, 1.0, 0.05)


#2, select an analysis method
#Home
if selected_option=="Home":
    image = Image.open('./Fig/0.png')
    st.image(image, caption='',use_column_width=True)
    st.subheader('1, Upload omics data')
    f=open('./Fig/1.txt', 'r')
    st.write(f.read()) 
    image = Image.open('./Fig/1.png')
    st.image(image, caption='',use_column_width=True)
    st.subheader('2, Select an analysis method')
    image = Image.open('./Fig/2.png')
    st.image(image, caption='',use_column_width=True)
    f=open('./Fig/2A.txt', 'r')
    st.write(f.read()) 
    f=open('./Fig/2B.txt', 'r')
    st.write(f.read()) 
    f=open('./Fig/2C.txt', 'r')
    st.write(f.read()) 
    f=open('./Fig/2D.txt', 'r')
    st.write(f.read()) 
    st.subheader('3, Set parameters (option)')
    f=open('./Fig/3.txt', 'r')
    st.write(f.read())
    f.close()
    st.subheader('License')
    f=open('./Fig/4.txt', 'r')
    st.write(f.read())
    f.close()
    
#Output
if selected_option=="A, gene regulatory network (including TF, miRNA, and mRNA) (transcriptome)":
    if Tran is not None:
        st.subheader('A, gene regulatory network (including TF, miRNA, and mRNA)')
        #Delete
        if(os.path.isfile('TFmiRNA-mRNA.zip')):
            os.remove('TFmiRNA-mRNA.zip')
        if(os.path.isfile('mRNA-mRNA.zip')):
            os.remove('mRNA-mRNA.zip')
        if(os.path.isfile('protein-protein.zip')):
            os.remove('protein-protein.zip') 
        if(os.path.isfile('mRNAmetabolite-Enzyme.zip')):
            os.remove('mRNAmetabolite-Enzyme.zip') 
        if(os.path.isfile('Transporter.zip')):
            os.remove('Transporter.zip')  
        #Database
        Name=pd.read_csv("./Database/230228Molecule2Name.csv")
        miRNA=pd.read_csv("./Database/230424mmu_miRNA_Gene.csv").drop_duplicates(subset=['Name', 'miRNA'])
        
        #mRNA_miRNA
        Tran2=pd.merge(Tran.rename(columns={'ENSMUSG': 'Molecule'}),Name, on='Molecule', how='inner')
        Tran3_UP=Tran2[Tran2["FC"]=="UP"].drop_duplicates(subset=['Name']).drop('Molecule', axis=1)
        Tran3_Down=Tran2[Tran2["FC"]=="Down"].drop_duplicates(subset=['Name']).drop('Molecule', axis=1)
        miRNAmRNA1=pd.concat([pd.merge(Tran3_UP, miRNA, on='Name', how='inner'),
                            pd.merge(Tran3_Down, miRNA, on='Name', how='inner')]).iloc[:,[2,1,0]]
        miRNAmRNA1["Regulation"]="miRNA"
        miRNAmRNA1.rename(columns={'miRNA': 'Regulator'}, inplace=True)
        miRNAmRNA1.rename(columns={'Name': 'mRNA'}, inplace=True)
        del miRNA
        del Name
        del Tran2
        gc.collect()
                
        #TF_mRNA        
        @st.cache_data
        def TFmiRNA_estimation():
            Scores=pd.concat([pd.read_csv("./Database/230228ChipScore1.csv", dtype="int16"),
                    pd.read_csv("./Database/230228ChipScore2.csv", dtype="int16"),
                    pd.read_csv("./Database/230228ChipScore3.csv", dtype="int16"),
                    pd.read_csv("./Database/230228ChipScore4.csv", dtype="int16"),
                    pd.read_csv("./Database/230228ChipScore5.csv", dtype="int16"),
                    pd.read_csv("./Database/230228ChipScore6.csv", dtype="int16")],axis=0).reset_index(drop=True)
            Scores=pd.concat([pd.read_csv("./Database/230228Genename.csv"),Scores],axis=1)   
            return Scores
        Scores=TFmiRNA_estimation()            
        TF=pd.read_csv("./Database/230228TF.csv")
      
        #TF_UP
        P = pd.DataFrame([], columns=['TF', 'p','list'], index=range(len(TF)))           
        for i in range (0,len(TF)):
            #DEG
            Cuery=pd.merge(Scores[['Name',str(i)]],Tran3_UP, how="inner", on = "Name",copy=False)
            Cuery_hit=Cuery[str(i)]>=Thre
            P["list"][i]=Cuery[Cuery[str(i)]>=Thre]["Name"].values
            #Background
            Cuery_ref=pd.merge(Scores[['Name',str(i)]],df_ref, how="inner", on = "Name",copy=False)
            Cuery_refhit=Cuery_ref[str(i)]>=Thre
            res = fisher_exact([[Cuery_hit.sum(),Cuery_refhit.sum()],
                                    [len(Tran3_UP)-Cuery_hit.sum(), len(df_ref)-Cuery_refhit.sum()]], 
                               alternative='greater')
            #stats.barnard_exact or stats.boschloo_exact are better, but take time
            P["TF"][i]=TF["TF"][i]
            P["p"][i]=res.pvalue
        P=pd.concat([P,pd.read_csv("./Database/230228experimentList_mm10.csv")[["Var4","Var5"]]],axis=1)
        if Organ=="All":
            P=P.reset_index(drop=True)
        else:
            P=P[P["Var5"]==Organ].reset_index(drop=True)
        Q=pd.DataFrame([multi.multipletests(P["p"], alpha=0.05, method='fdr_bh')[1]]).T
        P=pd.concat([P,Q],axis=1)
        Q3_UP=P[P[0]<FDR][['TF', 'p','list',0]].rename(columns={0: 'Q'})

        #TF_Down
        P = pd.DataFrame([], columns=['TF', 'p','list'], index=range(len(TF))) 
        for i in range (0,len(TF)):
            #DEG
            Cuery=pd.merge(Scores[['Name',str(i)]],Tran3_Down, how="inner", on = "Name",copy=False)
            Cuery_hit=Cuery[str(i)]>=Thre
            P["list"][i]=Cuery[Cuery[str(i)]>=Thre]["Name"].values
            #Background
            Cuery_ref=pd.merge(Scores[['Name',str(i)]],df_ref, how="inner", on = "Name",copy=False)
            Cuery_refhit=Cuery_ref[str(i)]>=Thre
            res = fisher_exact([[Cuery_hit.sum(),Cuery_refhit.sum()],
                                    [len(Tran3_Down)-Cuery_hit.sum(), len(df_ref)-Cuery_refhit.sum()]], alternative='greater')
            #stats.barnard_exact or stats.boschloo_exact are better, but take time
            P["TF"][i]=TF["TF"][i]
            P["p"][i]=res.pvalue
        P=pd.concat([P,pd.read_csv("./Database/230228experimentList_mm10.csv")[["Var4","Var5"]]],axis=1)
        if Organ=="All":
            P=P.reset_index(drop=True)
        else:
            P=P[P["Var5"]==Organ].reset_index(drop=True)
        Q=pd.DataFrame([multi.multipletests(P["p"], alpha=0.05, method='fdr_bh')[1]]).T
        P=pd.concat([P,Q],axis=1)
        Q3_Down=P[P[0]<FDR][['TF', 'p','list',0]].rename(columns={0: 'Q'})
        st.cache_data.clear()
        del Scores
        del TF
        gc.collect()
        
        #TF
        Q3_UP["Regulation"]="UP"
        Q3_Down["Regulation"]="Down"
        Q4=pd.concat([Q3_UP,Q3_Down]).reset_index(drop=True)
        del Q3_UP
        del Q3_Down
        gc.collect()        
        TFmRNA=pd.DataFrame([], columns=['TF','mRNA', 'FC'])
        for i in range (0,len(Q4)):
            TFmRNA1=pd.DataFrame([Q4["list"][i]]).T.rename(columns={0: 'mRNA'})
            TFmRNA1["TF"]=Q4["TF"][i]
            TFmRNA1["FC"]=Q4["Regulation"][i]
            TFmRNA=pd.concat([TFmRNA,TFmRNA1])
        TFmRNA["Regulation"]="TF"
        TFmRNA.rename(columns={'TF': 'Regulator'}, inplace=True)
        
        #TFmiRNA-mRNA
        TFmiRNAmRNA=pd.concat([TFmRNA,miRNAmRNA1])
        del Q4
        del TFmRNA
        del miRNAmRNA1
        gc.collect()  
        TFmiRNAmRNA=TFmiRNAmRNA.drop_duplicates(subset=['Regulator', 'mRNA','FC','Regulation'])
        TFmiRNAmRNAcopy=TFmiRNAmRNA.copy()
        TFmiRNAmRNAcopy.rename(columns={'Regulator': 'TF/miRNA'}, inplace=True)
        #st.subheader('TF/miRNA-mRNA')
        st.write(TFmiRNAmRNAcopy.set_index('TF/miRNA'))
        f=open('./Fig/A1.txt', 'r')
        st.write(f.read()) 
        #while not os.path.isfile('TFmiRNA-mRNA.zip'):
        if(os.path.isfile('TFmiRNA-mRNA.zip')):
            os.remove('TFmiRNA-mRNA.zip')
        with zipfile.ZipFile('TFmiRNA-mRNA.zip', 'x') as csv_zip:
            csv_zip.writestr("TFmiRNA-mRNA.csv",
                            TFmiRNAmRNAcopy.to_csv(index=False))
        with open("TFmiRNA-mRNA.zip", "rb") as file: 
            st.download_button(label = "Download TFmiRNA-mRNA data",data = file,file_name = "TFmiRNA-mRNA.zip")        
        
        del TFmiRNAmRNAcopy
        gc.collect()
        
        #Network
        fig = plt.figure(figsize=(12,8),facecolor='black')
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.set_facecolor('black')
        G1 = nx.Graph()
        G1.add_nodes_from(TFmiRNAmRNA["Regulator"].unique())
        G1.add_nodes_from(TFmiRNAmRNA["mRNA"].unique())
        edge_lists=[tuple(x) for x in TFmiRNAmRNA.iloc[:,[0,1]].values]
        G1.add_edges_from(edge_lists)

        pos = {}
        for i in G1.nodes():
            if (TFmiRNAmRNA[TFmiRNAmRNA["Regulator"]==i]["Regulation"]=="TF").sum()>0:
                random.seed(1) 
                pos[i] = (1.5*np.random.rand()-0.3,1.5*np.random.rand(),2) 
            elif (TFmiRNAmRNA[TFmiRNAmRNA["Regulator"]==i]["Regulation"]=="miRNA").sum()>0:
                random.seed(2) 
                pos[i] = (-1.5*np.random.rand()+0.3,-1.5*np.random.rand()-0.2,1)
            else:
                random.seed(3)
                pos[i] = (3*np.random.rand()-1.5,3*np.random.rand()-1.5,0) 
        nx.set_node_attributes(G1,pos,'pos')
        node_pos = nx.get_node_attributes(G1,'pos')
        x,y,z = zip(*node_pos.values())
        ax1.scatter(x,y,z,marker='o',s=0.001,c='black',alpha=0.5)

        for e in G1.edges():
            TFmiRNAmRNA1=TFmiRNAmRNA[TFmiRNAmRNA["Regulator"]==e[0]]
            TFmiRNAmRNA1=TFmiRNAmRNA1[TFmiRNAmRNA1["mRNA"]==e[1]]
            if (TFmiRNAmRNA1["FC"]=="UP").sum()>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.2,alpha=0.5)
            if (TFmiRNAmRNA1["FC"]=="Down").sum()>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.2,alpha=0.5)        

        verts = [[(-2.0, -2.0, 0), (-2.0, 2.0, 0), (2.0, 2.0, 0), (2.0, -2.0, 0)],
                [(-2.0, -2.0, 1), (-2.0, 0, 1), (2.0, 0, 1), (2.0, -2.0, 1)],
                [(-2.0, 0, 2), (-2.0, 2.0, 2), (2.0, 2.0, 2), (2.0, 0, 2)]]
        planes = Poly3DCollection(verts, alpha=0.3, facecolor='lightgray')

        ax1.text(-2.5,2.2,1.8,
                "TF ("+"I:"+str(len(TFmiRNAmRNA[(TFmiRNAmRNA["Regulation"]=="TF") & (TFmiRNAmRNA["FC"]=="UP")]["Regulator"].unique()))+', D:'+str(len(TFmiRNAmRNA[(TFmiRNAmRNA["Regulation"]=="TF") & (TFmiRNAmRNA["FC"]=="Down")]["Regulator"].unique()))+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,0.8,
                "miRNA ("+"I:"+str(len(TFmiRNAmRNA[(TFmiRNAmRNA["Regulation"]=="miRNA") & (TFmiRNAmRNA["FC"]=="UP")]["Regulator"].unique()))+', D:'+str(len(TFmiRNAmRNA[(TFmiRNAmRNA["Regulation"]=="miRNA") & (TFmiRNAmRNA["FC"]=="Down")]["Regulator"].unique()))+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,-0.2,
                "mRNA ("+"I:"+str(len(TFmiRNAmRNA[TFmiRNAmRNA['FC']=="UP"]["mRNA"].unique()))+', D:'+str(len(TFmiRNAmRNA[TFmiRNAmRNA['FC']=="Down"]["mRNA"].unique()))+')',
                fontname="Arial",color="white",fontsize=9)

        ax1.add_collection3d(planes)
        ax1.set_xlim([-1,4])
        ax1.set_ylim([-1,4])
        ax1.set_zlim([-1,4])
        ax1.axis('off')
        ax1.view_init(7, 10)
        st.pyplot(fig)
        f=open('./Fig/A2.txt', 'r')
        st.write(f.read()) 
        
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G1)
        for node in net.get_nodes():
            random.seed(4) 
            net.get_node(node)['x']=5000*pos[node][0]
            net.get_node(node)['y']=500*np.random.rand()-5000*pos[node][2]
            net.get_node(node)['physics']=False
            net.get_node(node)['label']=str(node) 
            net.get_node(node)['size']=5
            if len(TFmiRNAmRNA[TFmiRNAmRNA["Regulator"]==node])>0:
                net.get_node(node)['color']="white"
            if len(TFmiRNAmRNA[TFmiRNAmRNA["mRNA"]==node])>0:
                TFmiRNAmRNA1=TFmiRNAmRNA[TFmiRNAmRNA["mRNA"]==node]
                if (TFmiRNAmRNA1["FC"]=="UP").sum()>0:
                    net.get_node(node)['color']=Inc_color
                elif (TFmiRNAmRNA1["FC"]=="Down").sum()>0:
                    net.get_node(node)['color']=Dec_color
        for edge in net.get_edges():
            TFmiRNAmRNA1=TFmiRNAmRNA[TFmiRNAmRNA["Regulator"]==edge["from"]]
            TFmiRNAmRNA1=TFmiRNAmRNA1[TFmiRNAmRNA1["mRNA"]==edge['to']]
            if (TFmiRNAmRNA1["FC"]=="UP").sum()>0:
                edge['color']=Inc_color
            elif (TFmiRNAmRNA1["FC"]=="Down").sum()>0:
                edge['color']=Dec_color
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("TFmiRNA-mRNA.html")
        HtmlFile = open("TFmiRNA-mRNA.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("TFmiRNA-mRNA.html", 'r'),
                        file_name="TFmiRNA-mRNA.html")
        G1 = nx.Graph()
        del TFmiRNAmRNA
        del TFmiRNAmRNA1
        gc.collect()  
        
        current_variables = list(globals().keys())
        exclude_list = ['current_variables', 'exclude_list','selected_option']
        variables_to_delete = [var for var in current_variables if var not in exclude_list]

        for var_name in variables_to_delete:
            del globals()[var_name]
        import gc
        gc.collect()            
        
    else:
        st.write("Please upload transciptome data (organ or cell).")


if selected_option=="B, mRNA (protein)-mRNA (protein) interaction (transcriptome or proteome)":
    if Tran is not None:
        st.subheader('B, mRNA-mRNA interaction')
        #Delete
        if(os.path.isfile('TFmiRNA-mRNA.zip')):
            os.remove('TFmiRNA-mRNA.zip')
        if(os.path.isfile('mRNA-mRNA.zip')):
            os.remove('mRNA-mRNA.zip')
        if(os.path.isfile('protein-protein.zip')):
            os.remove('protein-protein.zip') 
        if(os.path.isfile('mRNAmetabolite-Enzyme.zip')):
            os.remove('mRNAmetabolite-Enzyme.zip') 
        if(os.path.isfile('Transporter.zip')):
            os.remove('Transporter.zip')  
        #Database
        Name=pd.read_csv("./Database/230228Molecule2Name.csv")
        PPI=pd.read_csv("./Database/230228String400_2.csv")
        Gene2Protein=pd.read_csv("./Database/230228Gene2Protein.csv")
        #Output
        Tran1=pd.merge(Tran, Name.rename(columns={'Molecule': 'ENSMUSG'}), on='ENSMUSG', how='inner',copy=False)
        DEP=pd.merge(Tran1, Gene2Protein, on='ENSMUSG', how='inner').drop('ENSMUSG', axis=1).drop_duplicates(subset=['FC', 'Name', 'ENSMUSP'])
        DEP3=pd.merge(DEP.rename(columns={'ENSMUSP': 'protein1'}), PPI, on='protein1', how='inner',copy=False)
        DEP3=pd.merge(DEP3, DEP.rename(columns={'ENSMUSP': 'protein2'}), on='protein2', how='inner',copy=False)
        DEP_UP=DEP3[(DEP3['FC_x'] =="UP") & (DEP3['FC_y'] =="UP")]
        DEP_Down=DEP3[(DEP3['FC_x'] =="Down") & (DEP3['FC_y'] =="Down")]
        del DEP
        del DEP3
        del Name
        gc.collect()
                
        #Network_Dat
        Gd = nx.Graph()
        Gd.add_nodes_from(pd.concat([PPI["protein1"],PPI["protein2"]]).unique())
        edge_lists=[tuple(x) for x in PPI[["protein1","protein2"]].values]
        Gd.add_edges_from(edge_lists)
        degree1 = list(dict(nx.degree(Gd)).values())
        x1 = [i for i in range(max(degree1)+1)]
        y1 = [degree1.count(i) for i in x1]
        Pk10=pd.DataFrame([x1,y1]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk11=Pk10[Pk10["k"]>0]
        Pk12=Pk11[Pk11["P(k)"]>0]
        del PPI
        gc.collect()
        
        #Network_UP
        G = nx.Graph()
        G.add_nodes_from(DEP_UP["Name_x"].unique(), bipartite=0)
        G.add_nodes_from(DEP_UP["Name_y"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in DEP_UP[["Name_x","Name_y"]].values]
        G.add_edges_from(edge_lists)
        degree = list(dict(nx.degree(G)).values())
        x = [i for i in range(max(degree)+1)]
        y = [degree.count(i) for i in x]
        Pk=pd.DataFrame([x,y]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk1=Pk[Pk["k"]>0]
        Pk2=Pk1[Pk1["P(k)"]>0]
        
        Centrality=pd.DataFrame([nx.degree_centrality(G)]).T.rename(columns={0: 'Degree centrality'})
        Centrality1=pd.DataFrame(index=Centrality.index,columns=['Degree centrality',
                                                                'Degree', 'Neighbors'])
        for i in range (len(Centrality)):
            Centrality1["Degree centrality"][i]=Centrality["Degree centrality"][i]
            Centrality1["Degree"][i]=G.degree(Centrality.index[i])
            Centrality1["Neighbors"][i]=list(nx.all_neighbors(G,Centrality.index[i]))
        st.write(Centrality1.sort_values('Degree', ascending=False))
        Increased_output=Centrality1.copy()
        Increased_output=Increased_output.sort_values('Degree', ascending=False)
        f=open('./Fig/B1.txt', 'r')
        st.write(f.read()) 
        
        st.write(pd.DataFrame(data=np.array([[nx.density(Gd),nx.average_clustering(Gd),nx.degree_pearson_correlation_coefficient(Gd,0),np.mean(degree1)],
                        [nx.density(G),nx.average_clustering(G),nx.degree_pearson_correlation_coefficient(G,0),np.mean(degree)]]),
                    index=['Database', 'Uploaded data'],columns=['Density', 'Clustering coefficient', 'Assortativity','Mean degree']))
        f=open('./Fig/B7.txt', 'r')
        st.write(f.read())
         
        ##degree distribution
        fig, ax = plt.subplots(figsize=(2, 2))
        #Database
        df_X =np.log10(Pk12["k"])
        df_y =np.log10(Pk12["P(k)"])
        color="gray"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Database:"+str(-(result.params["k"]))[0:4]+" ("+str(-(result.conf_int(alpha=0.05).T["k"][1]))[0:4]+" - "+str(-(result.conf_int(alpha=0.05).T["k"][0]))[0:4]+")")


        #Upload
        df_X =np.log10(Pk2["k"])
        df_y =np.log10(Pk2["P(k)"])
        color="black"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Upload:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel("log"+"${_1}$"+"${_0}$"+"k")
        ax.set_ylabel("log"+"${_1}$"+"${_0}$"+"N(k)")
        ax.set_title("Increased")
        ax.legend(bbox_to_anchor=(2.5, 0.5), loc="upper right",fontsize=8)
        st.pyplot(fig)
        f=open('./Fig/B2.txt', 'r')
        st.write(f.read()) 

        fig, ax = plt.subplots(figsize=(10,10),facecolor='black')
        ax.set_facecolor('black')
        pos = nx.spring_layout(G, k=0.1)
        #pos =nx.circular_layout(G)
        nx.draw_networkx_edges(G, pos, edge_color=Inc_color,width=0.5,ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=Inc_color,  node_size=2,alpha=1,ax=ax)
        ax.set_title("Increased (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")",color="white",fontsize=20)
        plt.axis('off')
        #st.pyplot(fig)

        st.subheader("Increased network (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")")
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G)
        for node in net.get_nodes():
            net.get_node(node)['x']=5000*pos[node][0]
            net.get_node(node)['y']=-5000*pos[node][1] #the minus is needed here to respect networkx y-axis convention 
            net.get_node(node)['physics']=False
            net.get_node(node)['label']=str(node) #set the node label as a string so that it can be displayed
            net.get_node(node)['color']=Inc_color
            net.get_node(node)['size']=5
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("Increased_PPI.html")
        HtmlFile = open("Increased_PPI.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("Increased_PPI.html", 'r'),
                        file_name="Increased_PPI.html")   
        f=open('./Fig/B3.txt', 'r')
        st.write(f.read())    
        
        #Network_Down
        G = nx.Graph()
        G.add_nodes_from(DEP_Down["Name_x"].unique(), bipartite=0)
        G.add_nodes_from(DEP_Down["Name_y"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in DEP_Down[["Name_x","Name_y"]].values]
        G.add_edges_from(edge_lists)
        degree = list(dict(nx.degree(G)).values())
        x = [i for i in range(max(degree)+1)]
        y = [degree.count(i) for i in x]
        Pk=pd.DataFrame([x,y]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk1=Pk[Pk["k"]>0]
        Pk2=Pk1[Pk1["P(k)"]>0]
        
        Centrality=pd.DataFrame([nx.degree_centrality(G)]).T.rename(columns={0: 'Degree centrality'})
        Centrality1=pd.DataFrame(index=Centrality.index,columns=['Degree centrality',
                                                                'Degree', 'Neighbors'])
        for i in range (len(Centrality)):
            Centrality1["Degree centrality"][i]=Centrality["Degree centrality"][i]
            Centrality1["Degree"][i]=G.degree(Centrality.index[i])
            Centrality1["Neighbors"][i]=list(nx.all_neighbors(G,Centrality.index[i]))
        st.write(Centrality1.sort_values('Degree', ascending=False))
        Decreased_output=Centrality1.copy()
        Decreased_output=Decreased_output.sort_values('Degree', ascending=False)

        f=open('./Fig/B4.txt', 'r')
        st.write(f.read())    

        st.write(pd.DataFrame(data=np.array([[nx.density(Gd),nx.average_clustering(Gd),nx.degree_pearson_correlation_coefficient(Gd,0),np.mean(degree1)],
                        [nx.density(G),nx.average_clustering(G),nx.degree_pearson_correlation_coefficient(G,0),np.mean(degree)]]),
                    index=['Database', 'Uploaded data'],columns=['Density', 'Clustering coefficient', 'Assortativity','Mean degree']))
        f=open('./Fig/B8.txt', 'r')
        st.write(f.read())

        ##degree distribution
        fig, ax = plt.subplots(figsize=(2, 2))
        #Database
        df_X =np.log10(Pk12["k"])
        df_y =np.log10(Pk12["P(k)"])
        color="gray"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Database:"+str(-(result.params["k"]))[0:4]+" ("+str(-(result.conf_int(alpha=0.05).T["k"][1]))[0:4]+" - "+str(-(result.conf_int(alpha=0.05).T["k"][0]))[0:4]+")")
                
       #Upload
        df_X =np.log10(Pk2["k"])
        df_y =np.log10(Pk2["P(k)"])
        color="black"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Upload:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel("log"+"${_1}$"+"${_0}$"+"k")
        ax.set_ylabel("log"+"${_1}$"+"${_0}$"+"N(k)")
        ax.set_title("Decreased")
        ax.legend(bbox_to_anchor=(2.5, 0.5), loc="upper right",fontsize=8)
        st.pyplot(fig)
        f=open('./Fig/B5.txt', 'r')
        st.write(f.read())   
        
        fig, ax = plt.subplots(figsize=(10,10),facecolor='black')
        ax.set_facecolor('black')
        pos = nx.spring_layout(G, k=0.1)
        #pos =nx.circular_layout(G)
        nx.draw_networkx_edges(G, pos, edge_color=Dec_color,width=0.5,ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=Dec_color,  node_size=2,alpha=1,ax=ax)
        ax.set_title("Decreased (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")",color="white",fontsize=20)
        plt.axis('off')
        #st.pyplot(fig)
        
        st.subheader("Decreased network (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")")
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G)
        for node in net.get_nodes():
            net.get_node(node)['x']=5000*pos[node][0]
            net.get_node(node)['y']=-5000*pos[node][1] 
            net.get_node(node)['physics']=False
            net.get_node(node)['label']=str(node)
            net.get_node(node)['color']=Dec_color
            net.get_node(node)['size']=5
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("Decreased_PPI.html")
        HtmlFile = open("Decreased_PPI.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("Decreased_PPI.html", 'r'),
                        file_name="Decreased_PPI.html")      
        Gd = nx.Graph()
        f=open('./Fig/B6.txt', 'r')
        st.write(f.read())
        if(os.path.isfile('mRNA-mRNA.zip')):
            os.remove('mRNA-mRNA.zip')   
        #while not os.path.isfile('mRNA-mRNA.zip'):
        with zipfile.ZipFile('mRNA-mRNA.zip', 'x') as csv_zip:
            csv_zip.writestr("Increased_mRNA-mRNA.csv",
                            Increased_output.to_csv(index=True))
            csv_zip.writestr("Decreased_mRNA-mRNA.csv",
                            Decreased_output.to_csv(index=True))                
        with open("mRNA-mRNA.zip", "rb") as file: 
            st.download_button(label = "Download mRNA-mRNA data",data = file,file_name = "mRNA-mRNA.zip")
    if Pro is not None:
        st.subheader('B, Protein-protein interaction')
        #Database
        Name=pd.read_csv("./Database/230228Molecule2Name.csv")
        PPI=pd.read_csv("./Database/230228String400_2.csv")
        #Gene2Protein=pd.read_csv("./Database/230228Gene2Protein.csv")
        #Output
        DEP=pd.merge(Pro, Name.rename(columns={'Molecule': 'ENSMUSP'}), on='ENSMUSP', how='inner',copy=False)
        DEP3=pd.merge(DEP.rename(columns={'ENSMUSP': 'protein1'}), PPI, on='protein1', how='inner',copy=False)
        DEP3=pd.merge(DEP3, DEP.rename(columns={'ENSMUSP': 'protein2'}), on='protein2', how='inner',copy=False)
        DEP_UP=DEP3[(DEP3['FC_x'] =="UP") & (DEP3['FC_y'] =="UP")]
        DEP_Down=DEP3[(DEP3['FC_x'] =="Down") & (DEP3['FC_y'] =="Down")]

        #Network_Dat
        Gd = nx.Graph()
        Gd.add_nodes_from(pd.concat([PPI["protein1"],PPI["protein2"]]).unique())
        edge_lists=[tuple(x) for x in PPI[["protein1","protein2"]].values]
        Gd.add_edges_from(edge_lists)
        degree1 = list(dict(nx.degree(Gd)).values())
        x1 = [i for i in range(max(degree1)+1)]
        y1 = [degree1.count(i) for i in x1]
        Pk10=pd.DataFrame([x1,y1]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk11=Pk10[Pk10["k"]>0]
        Pk12=Pk11[Pk11["P(k)"]>0]
        del Name
        del DEP3
        del PPI
        gc.collect()
        #Network_UP
        G = nx.Graph()
        G.add_nodes_from(DEP_UP["Name_x"].unique(), bipartite=0)
        G.add_nodes_from(DEP_UP["Name_y"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in DEP_UP[["Name_x","Name_y"]].values]
        G.add_edges_from(edge_lists)
        degree = list(dict(nx.degree(G)).values())
        x = [i for i in range(max(degree)+1)]
        y = [degree.count(i) for i in x]
        Pk=pd.DataFrame([x,y]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk1=Pk[Pk["k"]>0]
        Pk2=Pk1[Pk1["P(k)"]>0]
        
        Centrality=pd.DataFrame([nx.degree_centrality(G)]).T.rename(columns={0: 'Degree centrality'})
        Centrality1=pd.DataFrame(index=Centrality.index,columns=['Degree centrality',
                                                                'Degree', 'Neighbors'])
        for i in range (len(Centrality)):
            Centrality1["Degree centrality"][i]=Centrality["Degree centrality"][i]
            Centrality1["Degree"][i]=G.degree(Centrality.index[i])
            Centrality1["Neighbors"][i]=list(nx.all_neighbors(G,Centrality.index[i]))
        st.write(Centrality1.sort_values('Degree', ascending=False))
        Increased_output=Centrality1.copy()
        Increased_output=Increased_output.sort_values('Degree', ascending=False)
        f=open('./Fig/B1p.txt', 'r')
        st.write(f.read()) 
        

        st.write(pd.DataFrame(data=np.array([[nx.density(Gd),nx.average_clustering(Gd),nx.degree_pearson_correlation_coefficient(Gd,0),np.mean(degree1)],
                        [nx.density(G),nx.average_clustering(G),nx.degree_pearson_correlation_coefficient(G,0),np.mean(degree)]]),
                    index=['Database', 'Uploaded data'],columns=['Density', 'Clustering coefficient', 'Assortativity','Mean degree']))
        f=open('./Fig/B7.txt', 'r')
        st.write(f.read())
         
        ##degree distribution
        fig, ax = plt.subplots(figsize=(2, 2))
        #Database
        df_X =np.log10(Pk12["k"])
        df_y =np.log10(Pk12["P(k)"])
        color="gray"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Database:"+str(-(result.params["k"]))[0:4]+" ("+str(-(result.conf_int(alpha=0.05).T["k"][1]))[0:4]+" - "+str(-(result.conf_int(alpha=0.05).T["k"][0]))[0:4]+")")


        #Upload
        df_X =np.log10(Pk2["k"])
        df_y =np.log10(Pk2["P(k)"])
        color="black"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Upload:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel("log"+"${_1}$"+"${_0}$"+"k")
        ax.set_ylabel("log"+"${_1}$"+"${_0}$"+"N(k)")
        ax.set_title("Increased")
        ax.legend(bbox_to_anchor=(2.5, 0.5), loc="upper right",fontsize=8)
        st.pyplot(fig)
        f=open('./Fig/B2.txt', 'r')
        st.write(f.read()) 

        fig, ax = plt.subplots(figsize=(10,10),facecolor='black')
        ax.set_facecolor('black')
        pos = nx.spring_layout(G, k=0.1)
        #pos =nx.circular_layout(G)
        nx.draw_networkx_edges(G, pos, edge_color=Inc_color,width=0.5,ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=Inc_color,  node_size=2,alpha=1,ax=ax)
        ax.set_title("Increased (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")",color="white",fontsize=20)
        plt.axis('off')
        #st.pyplot(fig)

        st.subheader("Increased network (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")")
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G)
        for node in net.get_nodes():
            net.get_node(node)['x']=5000*pos[node][0]
            net.get_node(node)['y']=-5000*pos[node][1] #the minus is needed here to respect networkx y-axis convention 
            net.get_node(node)['physics']=False
            net.get_node(node)['label']=str(node) #set the node label as a string so that it can be displayed
            net.get_node(node)['color']=Inc_color
            net.get_node(node)['size']=5
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("Increased_PPI.html")
        HtmlFile = open("Increased_PPI.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("Increased_PPI.html", 'r'),
                        file_name="Increased_PPI.html")   
        f=open('./Fig/B3.txt', 'r')
        st.write(f.read())    
        
        #Network_Down
        G = nx.Graph()
        G.add_nodes_from(DEP_Down["Name_x"].unique(), bipartite=0)
        G.add_nodes_from(DEP_Down["Name_y"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in DEP_Down[["Name_x","Name_y"]].values]
        G.add_edges_from(edge_lists)
        degree = list(dict(nx.degree(G)).values())
        x = [i for i in range(max(degree)+1)]
        y = [degree.count(i) for i in x]
        Pk=pd.DataFrame([x,y]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk1=Pk[Pk["k"]>0]
        Pk2=Pk1[Pk1["P(k)"]>0]
        
        Centrality=pd.DataFrame([nx.degree_centrality(G)]).T.rename(columns={0: 'Degree centrality'})
        Centrality1=pd.DataFrame(index=Centrality.index,columns=['Degree centrality',
                                                                'Degree', 'Neighbors'])
        for i in range (len(Centrality)):
            Centrality1["Degree centrality"][i]=Centrality["Degree centrality"][i]
            Centrality1["Degree"][i]=G.degree(Centrality.index[i])
            Centrality1["Neighbors"][i]=list(nx.all_neighbors(G,Centrality.index[i]))
        st.write(Centrality1.sort_values('Degree', ascending=False))
        Decreased_output=Centrality1.copy()
        Decreased_output=Decreased_output.sort_values('Degree', ascending=False)

        f=open('./Fig/B4p.txt', 'r')
        st.write(f.read())    

        st.write(pd.DataFrame(data=np.array([[nx.density(Gd),nx.average_clustering(Gd),nx.degree_pearson_correlation_coefficient(Gd,0),np.mean(degree1)],
                        [nx.density(G),nx.average_clustering(G),nx.degree_pearson_correlation_coefficient(G,0),np.mean(degree)]]),
                    index=['Database', 'Uploaded data'],columns=['Density', 'Clustering coefficient', 'Assortativity','Mean degree']))
        f=open('./Fig/B8.txt', 'r')
        st.write(f.read())

        ##degree distribution
        fig, ax = plt.subplots(figsize=(2, 2))
        #Database
        df_X =np.log10(Pk12["k"])
        df_y =np.log10(Pk12["P(k)"])
        color="gray"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Database:"+str(-(result.params["k"]))[0:4]+" ("+str(-(result.conf_int(alpha=0.05).T["k"][1]))[0:4]+" - "+str(-(result.conf_int(alpha=0.05).T["k"][0]))[0:4]+")")
                
       #Upload
        df_X =np.log10(Pk2["k"])
        df_y =np.log10(Pk2["P(k)"])
        color="black"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Upload:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel("log"+"${_1}$"+"${_0}$"+"k")
        ax.set_ylabel("log"+"${_1}$"+"${_0}$"+"N(k)")
        ax.set_title("Decreased")
        ax.legend(bbox_to_anchor=(2.5, 0.5), loc="upper right",fontsize=8)
        st.pyplot(fig)
        f=open('./Fig/B5.txt', 'r')
        st.write(f.read())   
        
        fig, ax = plt.subplots(figsize=(10,10),facecolor='black')
        ax.set_facecolor('black')
        pos = nx.spring_layout(G, k=0.1)
        #pos =nx.circular_layout(G)
        nx.draw_networkx_edges(G, pos, edge_color=Dec_color,width=0.5,ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=Dec_color,  node_size=2,alpha=1,ax=ax)
        ax.set_title("Decreased (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")",color="white",fontsize=20)
        plt.axis('off')
        #st.pyplot(fig)
        
        st.subheader("Decreased network (Nodes:"+str(len(G.nodes))+", Edges:"+str(len(G.edges))+")")
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G)
        for node in net.get_nodes():
            net.get_node(node)['x']=5000*pos[node][0]
            net.get_node(node)['y']=-5000*pos[node][1] 
            net.get_node(node)['physics']=False
            net.get_node(node)['label']=str(node)
            net.get_node(node)['color']=Dec_color
            net.get_node(node)['size']=5
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("Decreased_PPI.html")
        HtmlFile = open("Decreased_PPI.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("Decreased_PPI.html", 'r'),
                        file_name="Decreased_PPI.html")      
        
        f=open('./Fig/B6.txt', 'r')
        st.write(f.read())
        Gd = nx.Graph()
        if(os.path.isfile('protein-protein.zip')):
            os.remove('protein-protein.zip')     
        #while not os.path.isfile('protein-protein.zip'):
        with zipfile.ZipFile('protein-protein.zip', 'x') as csv_zip:
            csv_zip.writestr("Increased_protein-protein.csv",
                            Increased_output.to_csv(index=True))
            csv_zip.writestr("Decreased_protein-protein.csv",
                            Decreased_output.to_csv(index=True))                
        with open("protein-protein.zip", "rb") as file: 
            st.download_button(label = "Download protein-protein data",data = file,file_name = "protein-protein.zip")              
        current_variables = list(globals().keys())
        exclude_list = ['current_variables', 'exclude_list','selected_option']
        variables_to_delete = [var for var in current_variables if var not in exclude_list]

        for var_name in variables_to_delete:
            del globals()[var_name]
        import gc
        gc.collect()    
    else:
        st.write("Please upload transciptome or proteome data (organ or cell).")


if selected_option=="C, metabolic network (including enzyme, mRNA, and metabolite) (transcriptome and metabolome (organ or cell))":
    if Meta1 is not None and Tran is not None:
        st.subheader('C, metabolic network (including enzyme, mRNA, and metabolite)')
        #Delete
        if(os.path.isfile('TFmiRNA-mRNA.zip')):
            os.remove('TFmiRNA-mRNA.zip')
        if(os.path.isfile('mRNA-mRNA.zip')):
            os.remove('mRNA-mRNA.zip')
        if(os.path.isfile('protein-protein.zip')):
            os.remove('protein-protein.zip') 
        if(os.path.isfile('mRNAmetabolite-Enzyme.zip')):
            os.remove('mRNAmetabolite-Enzyme.zip') 
        if(os.path.isfile('Transporter.zip')):
            os.remove('Transporter.zip')  
        #Database
        BRENDA= pd.read_csv("./Database/230228Metabolite2reaction_BRENDA.csv")
        Meta_KEGG= pd.read_csv("./Database/230228Metabolite2reaction_KEGG.csv")
        #Meta_KEGG2= pd.read_csv("./Database/230228Metabolite2reaction_KEGG.csv")
        Tran_KEGG= pd.read_csv("./Database/230228Transcript2reaction_KEGG.csv")
        Species=pd.read_csv("./Database/230228MMU.csv")
        Name=pd.read_csv("./Database/230228Molecule2Name.csv")
        MAP=pd.read_csv("./Database/230228MAP.csv")
        CPDA=pd.read_csv("./Database/MetabolomeTable_AdditionalID.csv")
        
        #CPD diplicate
        CPDA1=pd.merge(Meta1, CPDA, on='CPD', how='inner')[["CPD1","FC"]]
        CPDA1.rename(columns={'CPD1': 'CPD'}, inplace=True)
        Meta=pd.concat([Meta1,CPDA1]).reset_index(drop=True)

        #MMU only, metabolic pathway only
        BRENDA=pd.merge(Species[Species["MAP"]=="map01100"].rename(columns={"Molecule": 'Enzyme'}),
                BRENDA, on='Enzyme', how='inner',copy=False)
        Meta_KEGG=pd.merge(Species[Species["MAP"]=="map01100"].rename(columns={"Molecule": 'Enzyme'}),
                Meta_KEGG, on='Enzyme', how='inner',copy=False)
        Tran_KEGG=pd.merge(Species[Species["MAP"]=="map01100"].rename(columns={"Molecule": 'Enzyme'}),
                Tran_KEGG, on='Enzyme', how='inner',copy=False)

        #Connect (Metabolite/Transcript-Enzyme)
        BRENDA3=pd.merge(BRENDA, Meta, on='CPD', how='inner',copy=False).rename(columns={'CPD': 'Molecule'})
        BRENDA3["Regulation"]="Allosteric"
        Subpro=pd.merge(Meta_KEGG, Meta, on='CPD', how='inner',copy=False).rename(columns={'CPD': 'Molecule'})
        ALL=pd.concat([BRENDA3, Subpro])
        Tran_KEGG3=pd.merge(Tran_KEGG, Tran, on='ENSMUSG', how='inner',copy=False)
        Tran_KEGG3=Tran_KEGG3.drop("mmu", axis=1).rename(columns={'ENSMUSG': 'Molecule'})
        Tran_KEGG3["Regulation"]="Transcript"
        ALL=pd.concat([ALL, Tran_KEGG3])
        ALL=ALL.reset_index(drop=True)
        ALL=pd.merge(ALL, Name, on='Molecule', how='inner')
        Num=pd.DataFrame(ALL["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Total_number'})
        Num=Num.rename(columns={'index': 'Reaction'})

        #Metabolite
        A=pd.concat([BRENDA3, Subpro])
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Meta_number'})
        A_Num=A_Num.rename(columns={'index': 'Reaction'})
        Num=pd.merge(Num,A_Num,  on='Reaction', how='outer').fillna({'Meta_number': 0})
        #Transcript
        A_Num=pd.DataFrame(Tran_KEGG3["Enzyme"].value_counts()).reset_index().rename(columns={'index': 'Reaction'}).rename(columns={'Enzyme': 'Transcript_number'})
        Num=pd.merge(Num,A_Num,  on='Reaction', how='outer').fillna({'Transcript_number': 0})

        #Allosteric
        A=ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'UP') & (ALL['ACIN'] == 'AC')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Allo_AC_inc_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Allo_AC_inc_number': 0}).rename(columns={'Name': 'Allo_AC_inc_CPD'})

        A=ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'Down') & (ALL['ACIN'] == 'AC')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Allo_AC_dec_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Allo_AC_dec_number': 0}).rename(columns={'Name': 'Allo_AC_dec_CPD'})

        A=ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'UP') & (ALL['ACIN'] == 'IN')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Allo_IN_inc_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Allo_IN_inc_number': 0}).rename(columns={'Name': 'Allo_IN_inc_CPD'})

        A=ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'Down') & (ALL['ACIN'] == 'IN')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Allo_IN_dec_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Allo_IN_dec_number': 0}).rename(columns={'Name': 'Allo_IN_dec_CPD'})

        #Substrate
        A=ALL[(ALL['Regulation'] =="Substrate") & (ALL['FC'] == 'UP')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Sub_inc_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Sub_inc_number': 0}).rename(columns={'Name': 'Sub_inc_CPD'})

        A=ALL[(ALL['Regulation'] =="Substrate") & (ALL['FC'] == 'Down')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Sub_dec_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Sub_dec_number': 0}).rename(columns={'Name': 'Sub_dec_CPD'})

        #Product
        A=ALL[(ALL['Regulation'] =="Product") & (ALL['FC'] == 'UP')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Pro_inc_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Pro_inc_number': 0}).rename(columns={'Name': 'Pro_inc_CPD'})

        A=ALL[(ALL['Regulation'] =="Product") & (ALL['FC'] == 'Down')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Pro_dec_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Pro_dec_number': 0}).rename(columns={'Name': 'Pro_dec_CPD'})

        #Transcript
        A=ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'UP')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Tran_inc_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Tran_inc_number': 0}).rename(columns={'Name': 'Tran_inc_CPD'})

        A=ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'Down')]
        A1=A.groupby('Enzyme').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Enzyme"].value_counts()).reset_index().rename(columns={'Enzyme': 'Tran_dec_number'})
        AA=pd.merge(A_Num.rename(columns={'index': 'Reaction'}), 
                A1.rename(columns={'Enzyme': 'Reaction'}), on='Reaction', how='inner')
        Num=pd.merge(Num,AA,  on='Reaction', how='outer').fillna({'Tran_dec_number': 0}).rename(columns={'Name': 'Tran_dec_CPD'})

        #Pathway
        reaction=Species.groupby('Molecule').agg({'MAP Name': lambda x: ' | '.join(x)}).reset_index()
        Output_ALL1=pd.merge(reaction.rename(columns={'Molecule': 'Reaction'}),Num, 
                            on='Reaction', how='inner').sort_values('Total_number', ascending=False)
        st.write(Output_ALL1.set_index('Reaction'))
        output1=Output_ALL1.copy()
        f=open('./Fig/C1.txt', 'r')
        st.write(f.read())
         
        #Num. of Regulation
        Allo_dec_AC=len(ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'Down') & (ALL['ACIN'] == 'AC')])
        Allo_inc_AC=len(ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'UP') & (ALL['ACIN'] == 'AC')])
        Allo_dec_IN=len(ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'Down') & (ALL['ACIN'] == 'IN')])
        Allo_inc_IN=len(ALL[(ALL['Regulation'] =="Allosteric") & (ALL['FC'] == 'UP') & (ALL['ACIN'] == 'IN')])

        Substrate_act=len(ALL[(ALL['Regulation'] =="Substrate") & (ALL['FC'] == 'UP')])
        Substrate_inh=len(ALL[(ALL['Regulation'] =="Substrate") & (ALL['FC'] == 'Down')])
        Product_act=len(ALL[(ALL['Regulation'] =="Product") & (ALL['FC'] == 'UP')])
        Product_inh=len(ALL[(ALL['Regulation'] =="Product") & (ALL['FC'] == 'Down')])

        Trans_act=len(ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'UP')])
        Trans_inh=len(ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'Down')])

        Regulation=pd.DataFrame([[Trans_act, Trans_inh], [Substrate_act, Substrate_inh], 
                                [Product_act, Product_inh], [Allo_inc_AC, Allo_dec_AC], [Allo_inc_IN, Allo_dec_IN]],
                                columns=['Increased', 'decreased'],index=['Enzyme gene expression', 'Substrate', 'Product','Allosteric AC','Allosteric IN'])
        #display(Regulation)
        del A
        del A1
        del AA
        del Num
        gc.collect()        
        
        #Network
        fig = plt.figure(figsize=(12,8),facecolor='black')
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.set_facecolor('black')
        G1 = nx.Graph()
        G1.add_nodes_from(ALL[ALL["Regulation"]=="Transcript"]["Molecule"].unique())
        G1.add_nodes_from(ALL[ALL["Regulation"]!="Transcript"]["Molecule"].unique())
        G1.add_nodes_from(ALL["Enzyme"].unique())
        edge_lists=[tuple(x) for x in ALL[["Enzyme","Molecule"]].values]
        G1.add_edges_from(edge_lists)

        pos = {}
        for i in G1.nodes():
            if i[0]=="E":
                random.seed(1) 
                pos[i] = (3*np.random.rand()-1.5,3*np.random.rand()-1.5,2)
            elif i[0]=="C":
                random.seed(3) 
                pos[i] = (3*np.random.rand()-1.5,3*np.random.rand()-1.5,0)
            else:
                random.seed(2)
                pos[i] = (3*np.random.rand()-1.5,3*np.random.rand()-1.5,1)
        nx.set_node_attributes(G1,pos,'pos')
        node_pos = nx.get_node_attributes(G1,'pos')
        x,y,z = zip(*node_pos.values())
        ax1.scatter(x,y,z,marker='o',s=0.001,c='black',alpha=0.5)

        for e in G1.edges():
            ALL1=ALL[ALL["Enzyme"]==e[1]]
            ALL1=ALL1[ALL1["Molecule"]==e[0]]
            if len(ALL1[(ALL1["Regulation"]=="Transcript") & (ALL1['FC'] == 'UP')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Transcript") & (ALL1['FC'] == 'Down')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'UP') & (ALL1['ACIN'] == 'AC')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'Down') & (ALL1['ACIN'] == 'IN')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'Down') & (ALL1['ACIN'] == 'AC')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'UP') & (ALL1['ACIN'] == 'IN')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Substrate") & (ALL1['FC'] == 'UP')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Substrate") & (ALL1['FC'] == 'Down')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Product") & (ALL1['FC'] == 'Down')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1["Regulation"]=="Product") & (ALL1['FC'] == 'UP')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.4,alpha=0.5)
        verts = [list(zip([-0.5-1.5,-0.5-1.5,3.5-1.5,3.5-1.5],[-0.5-1.5,3.5-1.5,3.5-1.5,-0.5-1.5],[i,i,i,i])) for i in range(3)]
        planes = Poly3DCollection(verts, alpha=0.3, facecolor='lightgray')
        ax1.add_collection3d(planes)
        ax1.set_xlim([-1,4])
        ax1.set_ylim([-1,4])
        ax1.set_zlim([-1,4])
        ax1.axis('off')
        ax1.view_init(7, 10)

        ax1.text(-2.5,2.2,1.8,"Enzyme mRNA ("+"I:"+str(Regulation["Increased"][0])+', D:'+str(Regulation["decreased"][0])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,0.8,"Reaction ("+str(len(Output_ALL1))+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,0.4,"Substrate ("+"I:"+str(Regulation["Increased"][1])+', D:'+str(Regulation["decreased"][1])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,0.1,"Product ("+"I:"+str(Regulation["Increased"][2])+', D:'+str(Regulation["decreased"][2])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,-0.2,"Allosteric AC ("+"I:"+str(Regulation["Increased"][3])+', D:'+str(Regulation["decreased"][3])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,-0.5,"Allosteric IN ("+"I:"+str(Regulation["Increased"][4])+', D:'+str(Regulation["decreased"][4])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.2,-0.9,"Metabolite ("+"I:"+str(len(Meta1[Meta1["FC"]=="UP"]))+', D:'+str(len(Meta1[Meta1["FC"]=="Down"]))+')',
                fontname="Arial",color="white",fontsize=9)
        st.pyplot(fig)
        f=open('./Fig/C2.txt', 'r')
        st.write(f.read())
        
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G1)
        for node in net.get_nodes():
            random.seed(4) 
            net.get_node(node)['x']=5000*pos[node][1]
            net.get_node(node)['y']=500*np.random.rand()-5000*pos[node][2]
            net.get_node(node)['physics']=False
            net.get_node(node)['size']=5
            if len(ALL[ALL["Enzyme"]==node])>0:
                net.get_node(node)['label']=str(node)
                net.get_node(node)['color']="white"
            if len(ALL[ALL["Molecule"]==node])>0:
                net.get_node(node)['label']=Name[Name["Molecule"]==node]["Name"].unique()[0]
                ALL1=ALL[ALL["Molecule"]==node]
                if (ALL1["FC"]=="UP").sum()>0:
                    net.get_node(node)['color']=Inc_color
                elif (ALL1["FC"]=="Down").sum()>0:
                    net.get_node(node)['color']=Dec_color
        for edge in net.get_edges():
            ALL1=ALL[ALL["Enzyme"]==edge["to"]]
            ALL1=ALL1[ALL1["Molecule"]==edge["from"]]
            if len(ALL1[(ALL1["Regulation"]=="Transcript") & (ALL1['FC'] == 'UP')])>0:
                edge['color']=Inc_color
            elif len(ALL1[(ALL1["Regulation"]=="Transcript") & (ALL1['FC'] == 'Down')])>0:
                edge['color']=Dec_color
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'UP') & (ALL1['ACIN'] == 'AC')])>0:
                edge['color']=Inc_color
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'Down') & (ALL1['ACIN'] == 'IN')])>0:
                edge['color']=Inc_color
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'Down') & (ALL1['ACIN'] == 'AC')])>0:
                edge['color']=Dec_color
            elif len(ALL1[(ALL1["Regulation"]=="Allosteric") & (ALL1['FC'] == 'UP') & (ALL1['ACIN'] == 'IN')])>0:
                edge['color']=Dec_color
            elif len(ALL1[(ALL1["Regulation"]=="Substrate") & (ALL1['FC'] == 'UP')])>0:
                edge['color']=Inc_color
            elif len(ALL1[(ALL1["Regulation"]=="Substrate") & (ALL1['FC'] == 'Down')])>0:
                edge['color']=Dec_color
            elif len(ALL1[(ALL1["Regulation"]=="Product") & (ALL1['FC'] == 'Down')])>0:
                edge['color']=Inc_color
            elif len(ALL1[(ALL1["Regulation"]=="Product") & (ALL1['FC'] == 'UP')])>0:
                edge['color']=Dec_color
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("mRNAmetabolite-Enzyme.html")
        HtmlFile = open("mRNAmetabolite-Enzyme.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("mRNAmetabolite-Enzyme.html", 'r'),
                        file_name="mRNAmetabolite-Enzyme.html")           
        
        #Num. regulated reactions
        TraNum=len(Output_ALL1[Output_ALL1["Transcript_number"]>0])
        MetNum=len(Output_ALL1[Output_ALL1['Meta_number'] > 0])
        Both=TraNum+MetNum-len(Output_ALL1)
        Traonly=TraNum-Both
        Metaonly=MetNum-Both
        fig, ax = plt.subplots(figsize=(3,1.5))
        height = [Both, Traonly, Metaonly]  
        left = np.array([1, 2, 3])
        labels = ['Both'+" ("+str(Both)+")", 
                'Enzyme mRNA'+" ("+str(Traonly)+")", 'Metabolite'+" ("+str(Metaonly)+")"]
        ax.barh(left, height, color='gray', height=0.9, align='center')
        plt.yticks(left, labels)
        ax.set_xlabel("Num. of regulated reaction")
        ax.set_title("Metabolic pathway")
        st.pyplot(fig)
        f=open('./Fig/C3.txt', 'r')
        st.write(f.read())
            
        Path_all=pd.merge(Species.rename(columns={'Molecule': 'Enzyme'}),ALL.iloc[:,2:], on='Enzyme', how='inner')  
        #Percentage of regulations
        fig, ax = plt.subplots(figsize=(2, 2))
        for i in range(len(MAP["MAP Name"].unique())):
            Pathsub=Path_all[Path_all["MAP Name"]==MAP["MAP Name"].unique()[i]]
            TraNum=len(Pathsub[Pathsub["Regulation"]=='Transcript'].drop_duplicates(subset=['Enzyme']))
            MetNum=len(Pathsub[Pathsub["Regulation"]!='Transcript'].drop_duplicates(subset=['Enzyme']))
            TotalNum=len(Pathsub.drop_duplicates(subset=['Enzyme']))
            if TotalNum>0:
                ax.scatter(100*MetNum/TotalNum,100*TraNum/TotalNum, alpha=0.3, c="gray",s=TotalNum/10)
            if TotalNum>0.05*len(Path_all["Enzyme"].unique()):
                ax.text(100*MetNum/TotalNum,100*TraNum/TotalNum, MAP["MAP Name"].unique()[i], fontsize=3,color="black")
        ax.set_ylim([-10,110])
        ax.set_yticks([0,50,100])
        ax.set_xticks([0,50,100])
        ax.set_xlim([-10, 110])
        sns.despine()
        plt.xlabel("By Metabolite (%)")
        plt.ylabel("By Enzyme mRNA (%)")
        st.pyplot(fig)
        f=open('./Fig/C4.txt', 'r')
        st.write(f.read())
        G1 = nx.Graph()
        
        ##Network analysis
        Database_ALL=pd.concat([BRENDA.rename(columns={'CPD': 'Molecule'}),
                                Meta_KEGG.rename(columns={'CPD': 'Molecule'}),
                                Tran_KEGG.rename(columns={'ENSMUSG': 'Molecule'})])
        #Database_ALL
        G = nx.Graph()
        G.add_nodes_from(Database_ALL["Molecule"].unique(), bipartite=0)
        G.add_nodes_from(Database_ALL["Enzyme"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in Database_ALL[["Enzyme","Molecule"]].values]
        G.add_edges_from(edge_lists)
        Centrality=pd.DataFrame([nx.degree_centrality(G)]).T.rename(columns={0: 'Database_Degree centrality'})
        Centrality1=pd.DataFrame(index=Centrality.index,columns=['Database_Degree centrality',
                                                                'Database_Degree', 'Neighbors'])
        for i in range (len(Centrality)):
            Centrality1["Database_Degree centrality"][i]=Centrality["Database_Degree centrality"][i]
            Centrality1["Database_Degree"][i]=G.degree(Centrality.index[i])
            Centrality1["Neighbors"][i]=list(nx.all_neighbors(G,Centrality.index[i]))
        #st.write(Centrality1.sort_values('Database_Degree', ascending=False))
        #f=open('./Fig/C5.txt', 'r')
        #st.write(f.read())
        
        degree = list(dict(nx.degree(G)).values())
        x = [i for i in range(max(degree)+1)]
        y = [degree.count(i) for i in x]
        Pk=pd.DataFrame([x,y]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk1=Pk[Pk["k"]>0]
        Pk2=Pk1[Pk1["P(k)"]>0]
        del BRENDA
        del Meta_KEGG
        del Output_ALL1
        del Database_ALL
        gc.collect() 
        
        #Uploaded data
        G1 = nx.Graph()
        G1.add_nodes_from(ALL["Name"].unique(), bipartite=0)
        G1.add_nodes_from(ALL["Enzyme"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in ALL[["Enzyme","Name"]].values]
        G1.add_edges_from(edge_lists)
        Centrality=pd.DataFrame([nx.degree_centrality(G1)]).T.rename(columns={0: 'Uploaded_Degree centrality'})

        Centrality1=pd.DataFrame(index=Centrality.index,columns=['Uploaded_Degree centrality',
                                                                'Uploaded_Degree', 'Neighbors'])
        for i in range (len(Centrality)):
            Centrality1["Uploaded_Degree centrality"][i]=Centrality["Uploaded_Degree centrality"][i]
            Centrality1["Uploaded_Degree"][i]=G1.degree(Centrality.index[i])
            Centrality1["Neighbors"][i]=list(nx.all_neighbors(G1,Centrality.index[i]))
        st.write(Centrality1.sort_values('Uploaded_Degree', ascending=False))
        output2=Centrality1.copy()
        output2=output2.sort_values('Uploaded_Degree', ascending=False)
        f=open('./Fig/C6.txt', 'r')
        st.write(f.read())
        
        degree1 = list(dict(nx.degree(G1)).values())
        x1 = [i for i in range(max(degree1)+1)]
        y1 = [degree1.count(i) for i in x1]
        Pk10=pd.DataFrame([x1,y1]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk11=Pk10[Pk10["k"]>0]
        Pk12=Pk11[Pk11["P(k)"]>0]

        st.write(pd.DataFrame(data=np.array([[nx.density(G),nx.degree_pearson_correlation_coefficient(G,0),np.mean(degree)],
                        [nx.density(G1),nx.degree_pearson_correlation_coefficient(G1,0),np.mean(degree1)]]),
                    index=['Database', 'Uploaded data'],columns=['Density', 'Assortativity','Mean degree']))
        f=open('./Fig/C7.txt', 'r')
        st.write(f.read())
        
        #Database
        fig, ax = plt.subplots(figsize=(2, 2))
        df_X =np.log10(Pk2["k"])
        df_y =np.log10(Pk2["P(k)"])
        color="gray"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Database:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #Upload
        df_X =np.log10(Pk12["k"])
        df_y =np.log10(Pk12["P(k)"])
        color="black"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Upload:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel("log"+"${_1}$"+"${_0}$"+"k")
        ax.set_ylabel("log"+"${_1}$"+"${_0}$"+"N(k)")
        ax.legend(bbox_to_anchor=(2.5, 0.5), loc="upper right",fontsize=8)
        st.pyplot(fig)
        f=open('./Fig/C8.txt', 'r')
        st.write(f.read())
        G = nx.Graph()
        G1 = nx.Graph()
        #Loop
        G2 = nx.Graph()
        G2.add_nodes_from(ALL["Name"].unique(), bipartite=0)
        G2.add_nodes_from(ALL["Enzyme"].unique(), bipartite=1)
        edge_lists=[tuple(x) for x in ALL[["Enzyme","Name"]].values]
        G2.add_edges_from(edge_lists)
        ALL_L=pd.DataFrame(nx.cycle_basis(G2.to_undirected()))
        st.write(ALL_L)
        output3=ALL_L.copy()
        f=open('./Fig/C9.txt', 'r')
        st.write(f.read())
        del ALL_L
        gc.collect()                
        #Each pathway
        PathwayName = st.selectbox('Pathway:',MAP["MAP Name"].unique().tolist())
        i=np.where(MAP["MAP Name"].unique()==PathwayName)[0][0]
        #Make XML file
        #for i in range(0,len(MAP["MAP Name"].unique())):
        #    Pathsub=Path_all[Path_all["MAP Name"]==MAP["MAP Name"].unique()[i]].reset_index(drop=True)
        #    if len(Pathsub)>0:
        #        MAP1="mmu"+str(Pathsub["MAP"].unique()[0][3:])
        #        url = 'https://rest.kegg.jp/get/'+MAP1+'/kgml'
        #        response = urllib.request.urlopen(url)
        #        kgml = response.read().decode()
        #        root = ET.fromstring(kgml)
        #        with open('./xml/'+str(MAP["MAP Name"].unique()[i]).replace('/', '-')+'.xml', 'w') as f:
        #            f.write(ET.tostring(root, encoding='unicode'))
        Meta_KEGG2= pd.read_csv("./Database/230228Metabolite2reaction_KEGG.csv")
        Tran_KEGG8=Tran_KEGG.copy()
        Tran_KEGG8["mmu"]="mmu:"+ Tran_KEGG8['mmu'].astype(str)

        Pathsub=Path_all[Path_all["MAP Name"]==MAP["MAP Name"].unique()[i]].reset_index(drop=True)
        TraNum=len(Pathsub[Pathsub["Regulation"]=='Transcript'].drop_duplicates(subset=['Enzyme']))
        MetNum=len(Pathsub[Pathsub["Regulation"]!='Transcript'].drop_duplicates(subset=['Enzyme']))
        Both=TraNum+MetNum-len(Pathsub.drop_duplicates(subset=['Enzyme']))
        Traonly=TraNum-Both
        Metaonly=MetNum-Both

        fig, ax = plt.subplots(figsize=(2.5,1.5))
        height = [Both, Traonly, Metaonly]  
        left = np.array([1, 2, 3])
        labels = ['Both', 'Enzyme mRNA', 'Metabolite']
        ax.barh(left, height, color='gray', height=0.9, align='center')
        plt.yticks(left, labels)
        ax.set_xlabel("Num. of regulated reaction")
        ax.set_title(MAP["MAP Name"].unique()[i])
        st.pyplot(fig)
        f=open('./Fig/C10.txt', 'r')
        st.write(f.read())
                
        if len(Pathsub)>0:
            tree = ET.parse('./xml/'+str(MAP["MAP Name"].unique()[i]).replace('/', '-')+'.xml') 
            root = tree.getroot() 

            pos_dict = {}
            for entry in root.findall(".//entry"):
                graphics = entry.find('graphics')
                if str(type(graphics.get('x')))!= "<class 'NoneType'>":
                    node_id = entry.attrib['id']
                    node_name = entry.attrib['name']
                    node_type =entry.attrib['type']
                    x = float(graphics.get('x'))
                    y = -float(graphics.get('y'))
                    pos_dict[node_id] = (node_name,x, y)
            Path_loc=pd.DataFrame(pos_dict).T[0]
            Path_loc1=Path_loc.str.replace('cpd:', '')
            Path_loc2=Path_loc1.str.split(' ', expand=True)
            Path_loc3=pd.concat([pd.DataFrame(pos_dict).T[1].rename('x'),
                                    pd.DataFrame(pos_dict).T[2].rename('y'), Path_loc2],axis=1)

            #Enzyme (mmu>enzyme)
            m=0
            Enz_loc=pd.concat([Path_loc3["x"],Path_loc3["y"],pd.DataFrame([Path_loc3[m].rename('mmu')]).T],axis=1)
            Enz_loc1=pd.merge(Enz_loc, Tran_KEGG8, on='mmu', how='inner')
            for m in range (1,len(Path_loc3.columns)-2):
                Enz_loc=pd.concat([Path_loc3["x"],Path_loc3["y"],pd.DataFrame([Path_loc3[m].rename('mmu')]).T],axis=1)
                Enz_loc1=pd.concat([Enz_loc1,pd.merge(Enz_loc, Tran_KEGG8, on='mmu', how='inner')])
            Enz_loc2=Enz_loc1.drop(["mmu",'MAP',"MAP Name","ENSMUSG"], 
                                axis=1).drop_duplicates(subset=['Enzyme']).reset_index(drop=True)

            #DEM location
            Comp_loc=Path_loc3[Path_loc3[0].str.contains('C')].dropna(axis=1).rename(columns={0: 'Molecule'}).reset_index(drop=True)

            #Others
            Add=pd.merge(Comp_loc,Pathsub, on="Molecule", how ="outer", indicator=True).query(f'_merge == "right_only"').reset_index(drop=True)
            Add.loc[Add["Regulation"]=="Transcript", ['x']] = Path_loc3["x"].min()-100
            Add.loc[Add["Regulation"]!="Transcript", ['x']] = Path_loc3["x"].max()+100
            if len(Add)>0:
                #a=np.arange(-1200, -100,1000/len(Add))
                a=np.linspace(-1000, -100, len(Add))
                Addy=pd.DataFrame([a]).T.rename(columns={0: 'y'})
                Add=pd.concat([Add["Molecule"],Add["Enzyme"],Add["x"],Addy],axis=1).dropna()
                Add=Add.drop_duplicates(subset=['Enzyme','Molecule']).reset_index()

            G = nx.DiGraph()
            pos_dict = {}
            for j in range(len(Enz_loc2)):
                node_id=Enz_loc2["Enzyme"][j]
                pos_dict[node_id] = (Enz_loc2["x"][j], Enz_loc2["y"][j])
                G.add_node(node_id, label=Enz_loc2["Enzyme"][j])
            for j in range(len(Comp_loc)):
                node_id=Comp_loc["Molecule"][j]
                pos_dict[node_id] = (Comp_loc["x"][j], Comp_loc["y"][j])
                if len(Name[Name["Molecule"]==Comp_loc["Molecule"][j]]["Name"])>0:
                    node_name=Name[Name["Molecule"]==Comp_loc["Molecule"][j]]["Name"].str.split('; ', expand=True)[0].values[0]
                    node_name=node_name[:15]
                G.add_node(node_id, label=node_name)
            if len(Add)>0:
                for j in range(len(Add)):
                    node_id=Add["Molecule"][j]
                    pos_dict[node_id] = (Add["x"][j], Add["y"][j])
                    node_name=Name[Name["Molecule"]==Add["Molecule"][j]]["Name"].str.split('; ', expand=True)[0].values[0]
                    node_name=node_name[:15]
                    G.add_node(node_id, label=node_name)
                for j in range(len(Add)):
                    if len(Enz_loc2[Enz_loc2["Enzyme"]==Add["Enzyme"][j]])>0:
                        G.add_edge(Add["Molecule"][j],Add["Enzyme"][j])
            SubproPath=pd.merge(Species[Species["MAP Name"]==MAP["MAP Name"].unique()[i]].rename(columns={"Molecule": 'Enzyme'}),
                    Meta_KEGG2, on='Enzyme', how='inner')
            SubPro=pd.merge(SubproPath, Enz_loc2, on='Enzyme', how='inner')
            Node=pd.DataFrame([G.nodes]).T.reset_index().rename(columns={"index": 'CPD'})
            Sub=SubPro[SubPro["Regulation"]=="Substrate"].reset_index(drop=True)
            Sub1=pd.merge(Node, Sub, on='CPD', how='inner').drop_duplicates(subset=["CPD",'Enzyme']).reset_index(drop=True)
            Pro=SubPro[SubPro["Regulation"]=="Product"].reset_index(drop=True)
            Pro1=pd.merge(Node, Pro, on='CPD', how='inner').drop_duplicates(subset=["CPD",'Enzyme']).reset_index(drop=True)
            for j in range(len(Sub1)):
                G.add_edge(Sub1["CPD"][j],Sub1["Enzyme"][j])
            for j in range(len(Pro1)):
                G.add_edge(Pro1["Enzyme"][j],Pro1["CPD"][j])
            nodes_to_remove = [v for v in G if G.degree(v) < 1]
            G.remove_nodes_from(nodes_to_remove)
            #net=Network(notebook=True,height="750px",cdn_resources='remote' width="100%", bgcolor="#222222", font_color="white",directed=True)
            net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
                width="100%", bgcolor="#222222", font_color="white",directed=True)
            for node in G.nodes():
                x, y = pos_dict[node]
                size = 5
                if len(ALL[ALL["Enzyme"]==node])>0:
                    net.add_node(node, physics=False, x=x, y=-y, size=size,label=G.nodes[node]["label"],color="white",shape='diamond')
                if len(Enz_loc2[Enz_loc2["Enzyme"]==node])>0:
                    net.add_node(node, physics=False, x=x, y=-y, size=size,label=G.nodes[node]["label"],color="white",shape='diamond')
                elif len(ALL[ALL["Molecule"]==node])<1:
                    net.add_node(node, physics=False, x=x, y=-y, size=size,label=G.nodes[node]["label"],color="white")
                elif ALL[ALL["Molecule"]==node]["FC"].values[0]=="Down":
                    net.add_node(node, physics=False, x=x, y=-y, size=size,label=G.nodes[node]["label"],color=Dec_color,shape='triangleDown')
                elif ALL[ALL["Molecule"]==node]["FC"].values[0]=="UP":
                    net.add_node(node, physics=False, x=x, y=-y, size=size,label=G.nodes[node]["label"],color=Inc_color,shape='triangle')
            for e in G.edges():
                ALL1=Pathsub[Pathsub["Enzyme"]==e[0]]
                ALL1=ALL1[ALL1["Molecule"]==e[1]]
                ALL2=Pathsub[Pathsub["Enzyme"]==e[1]]
                ALL2=ALL2[ALL2["Molecule"]==e[0]]
                if len(ALL2[(ALL2["Regulation"]=="Transcript") & (ALL2['FC'] == 'UP')])>0:
                    edge_colors = Inc_color
                elif len(ALL2[(ALL2["Regulation"]=="Transcript") & (ALL2['FC'] == 'Down')])>0:
                    edge_colors = Dec_color
                elif len(ALL2[(ALL2["Regulation"]=="Allosteric") & (ALL2['FC'] == 'UP') & (ALL2['ACIN'] == 'AC')])>0:
                    edge_colors = Inc_color
                elif len(ALL2[(ALL2["Regulation"]=="Allosteric") & (ALL2['FC'] == 'Down') & (ALL2['ACIN'] == 'IN')])>0:
                    edge_colors = Inc_color
                elif len(ALL2[(ALL2["Regulation"]=="Allosteric") & (ALL2['FC'] == 'Down') & (ALL2['ACIN'] == 'AC')])>0:
                    edge_colors = Dec_color
                elif len(ALL2[(ALL2["Regulation"]=="Allosteric") & (ALL2['FC'] == 'UP') & (ALL2['ACIN'] == 'IN')])>0:
                    edge_colors =Dec_color
                elif len(ALL2[(ALL2["Regulation"]=="Substrate") & (ALL2['FC'] == 'UP')])>0:
                    edge_colors = Inc_color
                elif len(ALL2[(ALL2["Regulation"]=="Substrate") & (ALL2['FC'] == 'Down')])>0:
                    edge_colors = Dec_color
                elif len(ALL1[(ALL1["Regulation"]=="Product") & (ALL1['FC'] == 'Down')])>0:
                    edge_colors = Inc_color
                elif len(ALL1[(ALL1["Regulation"]=="Product") & (ALL1['FC'] == 'UP')])>0:
                    edge_colors = Dec_color
                else:
                    edge_colors = 'white'    
                net.add_edge(e[0], e[1], color=edge_colors,physics=False)
            #net.show_buttons() 
            net.show(str(MAP["MAP Name"].unique()[i]).replace('/', '-')+".html")
            HtmlFile = open(str(MAP["MAP Name"].unique()[i]).replace('/', '-')+".html", 'r')
            components.html(HtmlFile.read(), height=900)
            st.download_button(label="Download this network",data=open(str(MAP["MAP Name"].unique()[i]).replace('/', '-')+".html", 'r'),
                               file_name=str(MAP["MAP Name"].unique()[i]).replace('/', '-')+".html")   
            f=open('./Fig/C11.txt', 'r')
            st.write(f.read())  
            if(os.path.isfile('mRNAmetabolite-Enzyme.zip')):
                os.remove('mRNAmetabolite-Enzyme.zip')            
            with zipfile.ZipFile('mRNAmetabolite-Enzyme.zip', 'x') as csv_zip:
                csv_zip.writestr("mRNAmetabolite-Enzyme.csv",
                                output1.to_csv(index=False))
                csv_zip.writestr("Degree.csv",
                                output2.to_csv(index=True))
                csv_zip.writestr("List of loops.csv", 
                                output3.to_csv(index=False))
            with open("mRNAmetabolite-Enzyme.zip", "rb") as file: 
                st.download_button(label = "Download mRNAmetabolite-Enzyme data",data = file,file_name = "mRNAmetabolite-Enzyme.zip")
        del Meta_KEGG2
        #del Meta_KEGG
        del Tran_KEGG
        gc.collect()
        
        current_variables = list(globals().keys())
        exclude_list = ['current_variables', 'exclude_list','selected_option']
        variables_to_delete = [var for var in current_variables if var not in exclude_list]

        for var_name in variables_to_delete:
            del globals()[var_name]
        import gc
        gc.collect()    
    else:
        st.write("Please upload metabolome (organ or cell) and transciptome data (organ or cell).")
     
        
if selected_option=="D, metabolite exchange network (including transporter, mRNA, and metabolite) (transcriptome, metabolome (organ or cell) and metabolome (blood or medium))":
    if Meta1 is not None and Meta_blood is not None and Tran is not None:
        st.subheader('D, metabolic network (including transporter, mRNA, and metabolite)')
        #Delete
        if(os.path.isfile('TFmiRNA-mRNA.zip')):
            os.remove('TFmiRNA-mRNA.zip')
        if(os.path.isfile('mRNA-mRNA.zip')):
            os.remove('mRNA-mRNA.zip')
        if(os.path.isfile('protein-protein.zip')):
            os.remove('protein-protein.zip') 
        if(os.path.isfile('mRNAmetabolite-Enzyme.zip')):
            os.remove('mRNAmetabolite-Enzyme.zip') 
        if(os.path.isfile('Transporter.zip')):
            os.remove('Transporter.zip')
        
        #Input  
        Meta = Meta1.copy()
        Meta_Old = Meta1.copy()
        Meta_blood_Old=Meta_blood.copy()

        Name=pd.read_csv("./Database/230228Molecule2Name.csv")
        #CPD=pd.read_csv("./Database/230228CPD2Transporter.csv")
        ENSMUSG=pd.read_csv("./Database/230228ENSMUSG2Transporter.csv")
        CPDA=pd.read_csv("./Database/MetabolomeTable_AdditionalID.csv")
        Family=pd.read_csv("./Database/230228TransporterFamily.csv")

        #CPD duplicate
        CPDA1=pd.merge(Meta, CPDA, on='CPD', how='inner')
        CPDA1=CPDA1[["CPD1","FC"]]
        CPDA1.rename(columns={'CPD1': 'CPD'}, inplace=True)
        Meta=pd.concat([Meta,CPDA1]).reset_index(drop=True)

        CPDA1=pd.merge(Meta_blood, CPDA, on='CPD', how='inner')
        CPDA1=CPDA1[["CPD1","FC"]]
        CPDA1.rename(columns={'CPD1': 'CPD'}, inplace=True)
        Meta_blood=pd.concat([Meta_blood,CPDA1]).reset_index(drop=True)

        #Connect (Metabolite/Transcript-Transporter)
        OrganMeta=pd.merge(pd.read_csv("./Database/230228CPD2Transporter.csv"),Meta, on='CPD', how='inner',copy=False)
        OrganMeta.rename(columns={'CPD': 'Molecule'}, inplace=True)
        OrganMeta=pd.merge(OrganMeta,Name, on='Molecule', how='inner')
        OrganMeta["Molecule"]='Organ_' + OrganMeta['Molecule'].astype(str)
        OrganMeta["Regulation"]="Organ"
        BloodMeta=pd.merge(pd.read_csv("./Database/230228CPD2Transporter.csv"),Meta_blood, on='CPD', how='inner',copy=False)
        BloodMeta.rename(columns={'CPD': 'Molecule'}, inplace=True)
        BloodMeta=pd.merge(BloodMeta,Name, on='Molecule', how='inner',copy=False)
        BloodMeta["Molecule"]='Blood_' + BloodMeta['Molecule'].astype(str)
        BloodMeta["Regulation"]="Blood"
        OrganTran=pd.merge(ENSMUSG,Tran, on='ENSMUSG', how='inner')
        OrganTran.rename(columns={'ENSMUSG': 'Molecule'}, inplace=True)
        OrganTran=pd.merge(OrganTran,Name, on='Molecule', how='inner')
        OrganTran["Regulation"]="Transcript"
        ALL=pd.concat([OrganMeta, BloodMeta])
        ALL=pd.concat([ALL, OrganTran])
        ALL=ALL.reset_index(drop=True)
        Num=pd.DataFrame(ALL["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': 'Total_number'})
        Num=Num.rename(columns={'index': 'Transporter'})

        #Transporter-MetaboliteTranscript
        OrganMeta_Num=pd.DataFrame(OrganMeta["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': 'OrganMeta_number'})
        OrganMeta_Num=OrganMeta_Num.rename(columns={'index': 'Transporter'})
        Num=pd.merge(Num,OrganMeta_Num,  on='Transporter', how='outer').fillna({'OrganMeta_number': 0})

        BloodMeta_Num=pd.DataFrame(BloodMeta["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': 'BloodMeta_number'})
        BloodMeta_Num=BloodMeta_Num.rename(columns={'index': 'Transporter'})
        Num=pd.merge(Num,BloodMeta_Num,  on='Transporter', how='outer').fillna({'BloodMeta_number': 0})

        OrganTran_Num=pd.DataFrame(OrganTran["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': 'OrganTran_number'})
        OrganTran_Num=OrganTran_Num.rename(columns={'index': 'Transporter'})
        Num=pd.merge(Num,OrganTran_Num,  on='Transporter', how='outer').fillna({'OrganTran_number': 0})

        #Organ_Metabolite
        Regul="Organ"
        FC="UP"
        A=ALL[(ALL['FC'] == FC) & (ALL['Regulation'] == Regul)]
        A1=A.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': Regul+'_Metabo_'+FC+'_number'})
        A=pd.merge(A_Num.rename(columns={'index': 'Transporter'}), 
                A1, on='Transporter', how='inner')
        Num=pd.merge(Num,A,on='Transporter', how='outer').fillna({Regul+'_Metabo_'+FC+'_number': 0}).rename(columns={'Name': Regul+'_Metabo_'+FC+'_CPD'})
        Regul="Organ"
        FC="Down"
        A=ALL[(ALL['FC'] == FC) & (ALL['Regulation'] == Regul)]
        A1=A.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': Regul+'_Metabo_'+FC+'_number'})
        A=pd.merge(A_Num.rename(columns={'index': 'Transporter'}), 
                A1, on='Transporter', how='inner')
        Num=pd.merge(Num,A,on='Transporter', how='outer').fillna({Regul+'_Metabo_'+FC+'_number': 0}).rename(columns={'Name': Regul+'_Metabo_'+FC+'_CPD'})
        #Blood_Metabolite
        Regul="Blood"
        FC="UP"
        A=ALL[(ALL['FC'] == FC) & (ALL['Regulation'] == Regul)]
        A1=A.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': Regul+'_Metabo_'+FC+'_number'})
        A=pd.merge(A_Num.rename(columns={'index': 'Transporter'}), 
                A1, on='Transporter', how='inner')
        Num=pd.merge(Num,A,on='Transporter', how='outer').fillna({Regul+'_Metabo_'+FC+'_number': 0}).rename(columns={'Name': Regul+'_Metabo_'+FC+'_CPD'})
        Regul="Blood"
        FC="Down"
        A=ALL[(ALL['FC'] == FC) & (ALL['Regulation'] == Regul)]
        A1=A.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': Regul+'_Metabo_'+FC+'_number'})
        A=pd.merge(A_Num.rename(columns={'index': 'Transporter'}), 
                A1, on='Transporter', how='inner')
        Num=pd.merge(Num,A,on='Transporter', how='outer').fillna({Regul+'_Metabo_'+FC+'_number': 0}).rename(columns={'Name': Regul+'_Metabo_'+FC+'_CPD'})
        #Transcript
        Regul="Transcript"
        FC="UP"
        A=ALL[(ALL['FC'] == FC) & (ALL['Regulation'] == Regul)]
        A1=A.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': Regul+'_'+FC+'_number'})
        A=pd.merge(A_Num.rename(columns={'index': 'Transporter'}), 
                A1, on='Transporter', how='inner')
        Num=pd.merge(Num,A,on='Transporter', how='outer').fillna({Regul+'_'+FC+'_number': 0}).rename(columns={'Name': Regul+'_'+FC+'_ENSUMSG'})
        Regul="Transcript"
        FC="Down"
        A=ALL[(ALL['FC'] == FC) & (ALL['Regulation'] == Regul)]
        A1=A.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        A_Num=pd.DataFrame(A["Transporter"].value_counts()).reset_index().rename(columns={'Transporter': Regul+'_'+FC+'_number'})
        A=pd.merge(A_Num.rename(columns={'index': 'Transporter'}), 
                A1, on='Transporter', how='inner')
        Num=pd.merge(Num,A,on='Transporter', how='outer').fillna({Regul+'_'+FC+'_number': 0}).rename(columns={'Name': Regul+'_'+FC+'_ENSUMSG'})

        #Transporter name TransporterFamily
        Family1=pd.merge(Family,Num, on='Transporter', how='inner')
        Num=pd.merge(Family1[["Transporter","Transporter_Superfamily"]],Num, on='Transporter', how='outer')

        ENSMUSG.rename(columns={'ENSMUSG': 'Molecule'}, inplace=True)
        ENSMUSG1=pd.merge(ENSMUSG,Name, on='Molecule', how='inner')
        ENSMUSG2=ENSMUSG1.groupby('Transporter').agg({'Name': lambda x: ' | '.join(x)}).reset_index()
        ENSMUSG2.rename(columns={'Name': 'Transporter_Name'}, inplace=True)
        Num=pd.merge(ENSMUSG2,Num, on='Transporter', how='inner')
        Num=Num.sort_values('Total_number', ascending=False)
        st.write(Num.set_index('Transporter'))
        
        f=open('./Fig/D1.txt', 'r')
        st.write(f.read())    
        
        #Num. of Regulation
        Organ_Meta_UP=len(ALL[(ALL['Regulation'] =="Organ") & (ALL['FC'] == 'UP') ])
        Organ_Meta_Down=len(ALL[(ALL['Regulation'] =="Organ") & (ALL['FC'] == 'Down')])
        Blood_Meta_UP=len(ALL[(ALL['Regulation'] =="Blood") & (ALL['FC'] == 'UP') ])
        Blood_Meta_Down=len(ALL[(ALL['Regulation'] =="Blood") & (ALL['FC'] == 'Down')])
        Transcript_UP=len(ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'UP') ])
        Transcript_Down=len(ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'Down') ])

        Meta_Org_Inc_Num=len(ALL[(ALL['Regulation'] =="Organ") & (ALL['FC'] == 'UP')]["Molecule"].unique())
        Meta_Org_Dec_Num=len(ALL[(ALL['Regulation'] =="Organ") & (ALL['FC'] == 'Down')]["Molecule"].unique())
        Meta_blood_Inc_Num=len(ALL[(ALL['Regulation'] =="Blood") & (ALL['FC'] == 'UP')]["Molecule"].unique())
        Meta_blood_Dec_Num=len(ALL[(ALL['Regulation'] =="Blood") & (ALL['FC'] == 'Down')]["Molecule"].unique())
        Tran_Inc_Num=len(ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'UP')]["Molecule"].unique())
        Tran_Dec_Num=len(ALL[(ALL['Regulation'] =="Transcript") & (ALL['FC'] == 'Down')]["Molecule"].unique())

        Regulation=pd.DataFrame([[Tran_Inc_Num, Tran_Dec_Num], [Transcript_UP, Transcript_Down], [Organ_Meta_UP, Organ_Meta_Down], 
                                [Blood_Meta_UP, Blood_Meta_Down],[Meta_Org_Inc_Num, Meta_Org_Dec_Num],[Meta_blood_Inc_Num, Meta_blood_Dec_Num]],
                                columns=['Increased', 'decreased'],
                                index=['Enzyme gene expression','Enzyme gene expression_Regulation', 'Organ_Regulation', 'Blood_Regulation','Organ', 'Blood'])

        #Network
        fig = plt.figure(figsize=(12,8),facecolor='black')
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.set_facecolor('black')
        G1 = nx.Graph()
        G1.add_nodes_from(ALL["Molecule"].unique())
        G1.add_nodes_from(ALL["Transporter"].unique())
        edge_lists=[tuple(x) for x in ALL[["Molecule","Transporter"]].values]
        G1.add_edges_from(edge_lists)

        pos = {}
        for i in G1.nodes():
            if i[0]=="E":
                random.seed(1) 
                pos[i] = (3*np.random.rand()-1.5,3*np.random.rand()-1.5,2)
            elif i[0]=="O":
                random.seed(3) 
                pos[i] = (-1.5*np.random.rand()+1.0,-1.5*np.random.rand()-0.8,0)
            elif i[0]=="B":
                random.seed(3) 
                pos[i] = (1.5*np.random.rand()-1.0,1.5*np.random.rand()+0.8,-1)
            else:
                random.seed(2)
                pos[i] = (3*np.random.rand()-1.5,3*np.random.rand()-1.5,1)
        nx.set_node_attributes(G1,pos,'pos')
        node_pos = nx.get_node_attributes(G1,'pos')
        x,y,z = zip(*node_pos.values())
        ax1.scatter(x,y,z,marker='o',s=0.001,c='black',alpha=0.5)

        for e in G1.edges():
            ALL1=ALL[ALL["Transporter"]==e[1]]
            ALL1=ALL1[ALL1["Molecule"]==e[0]]
            if len(ALL1[(ALL1['FC'] == 'UP')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Inc_color,linewidth=0.4,alpha=0.5)
            elif len(ALL1[(ALL1['FC'] == 'Down')])>0:
                ax1.plot([node_pos[e[0]][0],node_pos[e[1]][0]],
                        [node_pos[e[0]][1],node_pos[e[1]][1]],
                        [node_pos[e[0]][2],node_pos[e[1]][2]],c=Dec_color,linewidth=0.4,alpha=0.5)
                
        verts = [[(-2.0, -2.0, 2), (-2.0, 2.0, 2), (2.0, 2.0, 2), (2.0, -2.0, 2)],
                [(-2.0, -2.0, 1), (-2.0, 2.0, 1), (2.0, 2.0, 1), (2.0, -2.0, 1)],
                [(-2.0, -2.5, 0), (-2.0, -0.5, 0), (2.0, -0.5, 0), (2.0, -2.5, 0)],
                [(-2.0, 0.5, -1), (-2.0, 2.5, -1), (2.0, 2.5, -1), (2.0, 0.5, -1)]]
        #verts = [[(-2.0, -2.0, 2), (-2.0, 2.0, 2), (2.0, 2.0, 2), (2.0, -2.0, 2)],
        #         [(-2.0, -2.0, 1), (-2.0, 2.0, 1), (2.0, 2.0, 1), (2.0, -2.0, 1)],
        #         [(-2.0, -2.0, 0), (-2.0, 0, 0), (2.0, 0, 0), (2.0, -2.0, 0)],
        #         [(-2.0, 0, -1), (-2.0, 2.0, -1), (2.0, 2, -1), (2.0, 0, -1)]]
        planes = Poly3DCollection(verts, alpha=0.3, facecolor='lightgray')
        ax1.add_collection3d(planes)
        ax1.set_xlim([-1,4])
        ax1.set_ylim([-1,4])
        ax1.set_zlim([-1,4])
        ax1.axis('off')
        ax1.view_init(7, 10)

        Meta_Old2=pd.merge(pd.read_csv("./Database/230228CPD2Transporter.csv"),Meta_Old, 
                           on='CPD', how='inner',copy=False)
        Meta_blood_Old2=pd.merge(pd.read_csv("./Database/230228CPD2Transporter.csv"),Meta_blood_Old,
                                 on='CPD', how='inner',copy=False)

        ax1.text(-2.5,2.0,1.8,"Transporter mRNA ("+"I:"+str(Regulation["Increased"][0])+', D:'+str(Regulation["decreased"][0])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.0,0.8,"Transporter ("+str(len(ALL["Transporter"].unique()))+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.0,0.3,"Organ_Transporter ("+"I:"+str(Regulation["Increased"][2])+', D:'+str(Regulation["decreased"][2])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.0,0,"Blood_Transporter ("+"I:"+str(Regulation["Increased"][3])+', D:'+str(Regulation["decreased"][3])+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.0,-0.5,
                "Metabolite (Organ) ("+"I:"+str(len(Meta_Old2[Meta_Old2["FC"]=="UP"]["CPD"].unique()))+', D:'+str(len(Meta_Old2[Meta_Old2["FC"]=="Down"]["CPD"].unique()))+')',
                fontname="Arial",color="white",fontsize=9)
        ax1.text(-2.5,2.0,-0.8,
                "Metabolite (Blood) ("+"I:"+str(len(Meta_blood_Old2[Meta_blood_Old2["FC"]=="UP"]["CPD"].unique()))+', D:'+str(len(Meta_blood_Old2[Meta_blood_Old2["FC"]=="Down"]["CPD"].unique()))+')',
                fontname="Arial",color="white",fontsize=9)
        st.pyplot(fig)
        f=open('./Fig/D2.txt', 'r')
        st.write(f.read())
        net=Network(notebook=True,cdn_resources='in_line',height="1000px", filter_menu=True, 
        width="100%", bgcolor="#222222", font_color="white",directed=False)
        net.from_nx(G1)
        for node in net.get_nodes():
            random.seed(4) 
            net.get_node(node)['x']=5000*pos[node][1]
            net.get_node(node)['y']=500*np.random.rand()-5000*pos[node][2] #the minus is needed here to respect networkx y-axis convention 
            net.get_node(node)['physics']=False
            net.get_node(node)['size']=5
            if len(ALL[ALL["Molecule"]==node])>0:
                net.get_node(node)['label']=ALL[ALL["Molecule"]==node]["Name"].unique()[0]
            else:
                net.get_node(node)['label']=str(node) 
            if pos[node][2]==2:
                ALL1=ALL[ALL["Molecule"]==node]
                if len(ALL1[ALL1["FC"]=="UP"])>0:
                    net.get_node(node)['color']=Inc_color
                else:
                    net.get_node(node)['color']=Dec_color
            if pos[node][2]==1:
                net.get_node(node)['color']="white"
            elif pos[node][2]==0:
                ALL1=ALL[ALL["Molecule"]==node]
                if len(ALL1[ALL1["FC"]=="UP"])>0:
                    net.get_node(node)['color']=Inc_color
                else:
                    net.get_node(node)['color']=Dec_color
            elif pos[node][2]==-1:
                ALL1=ALL[ALL["Molecule"]==node]
                if len(ALL1[ALL1["FC"]=="UP"])>0:
                    net.get_node(node)['color']=Inc_color
                else:
                    net.get_node(node)['color']=Dec_color   

        for edge in net.get_edges():
            edge['color']="white"
            ALL1=ALL[ALL["Molecule"]==edge["from"]]
            ALL1=ALL1[ALL1["Transporter"]==edge['to']]
            if (ALL1["FC"]=="UP").sum()>0:
                edge['color']=Inc_color
            elif (ALL1["FC"]=="Down").sum()>0:
                edge['color']=Dec_color
        net.toggle_physics(False)
        net.show_buttons() 
        net.show("Transporter.html")

        HtmlFile = open("Transporter.html", 'r')
        components.html(HtmlFile.read(), height=900)
        st.download_button(label="Download the interactice network",data=open("Transporter.html", 'r'),
                        file_name="Transporter.html") 
              
        fig, ax = plt.subplots(figsize=(3, 3))
        A=len(Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] == 0)])
        B=len(Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] > 0)])
        C=len(Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] > 0)])
        D=len(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] == 0)])
        E=len(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] == 0)])
        F=len(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] > 0)])
        G=len(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] > 0)])
        v=venn3(subsets = (A,B,C,D,E,F,G),
                set_labels = ('Metabolite (Blood)','Metabolite (Organ)', 'Enzyme mRNA'),
                set_colors=("gray", "gray", "gray"))
        st.pyplot(fig) 
        f=open('./Fig/D3.txt', 'r')
        st.write(f.read())
        

        st.write("Metabolite (Blood) & Metabolite (Organ) & Enzyme mRNA")
        st.write(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] > 0)].set_index('Transporter'))      
        st.write("Metabolite (Blood) & Metabolite (Organ) only")
        st.write(Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] > 0)].set_index('Transporter'))      
        st.write("Metabolite (Blood) & Enzyme mRNA only")
        st.write(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] == 0)].set_index('Transporter'))      
        st.write("Enzyme mRNA & Metabolite (Organ) only")
        st.write(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] > 0)].set_index('Transporter'))      
        st.write("Metabolite (Blood) only")
        st.write(Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] == 0)].set_index('Transporter'))      
        st.write("Enzyme mRNA only")
        st.write(Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] == 0)].set_index('Transporter'))      
        st.write("Metabolite (Organ) only")
        st.write(Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] > 0)].set_index('Transporter'))      
    
        #while not os.path.isfile('Transporter.zip'):
        if(os.path.isfile('Transporter.zip')):
            os.remove('Transporter.zip') 
        with zipfile.ZipFile('Transporter.zip', 'x') as csv_zip:
            csv_zip.writestr("Transporter.csv",
                            Num.to_csv(index=False))
            csv_zip.writestr("Metabolite (Blood) & Metabolite (Organ) & Enzyme mRNA.csv",
                            Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] > 0)].to_csv(index=False))
            csv_zip.writestr("Metabolite (Blood) & Metabolite (Organ) only.csv", 
                            Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] > 0)].to_csv(index=False))
            csv_zip.writestr("Metabolite (Blood) & Enzyme mRNA only.csv", 
                            Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] == 0)].to_csv(index=False)) 
            csv_zip.writestr("Enzyme mRNA & Metabolite (Organ) only.csv", 
                            Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] > 0)].to_csv(index=False))
            csv_zip.writestr("Metabolite (Blood) only.csv", 
                            Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] > 0)& (Num['OrganMeta_number'] == 0)].to_csv(index=False))        
            csv_zip.writestr("Enzyme mRNA only.csv", 
                            Num[(Num['OrganTran_number'] >0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] == 0)].to_csv(index=False))
            csv_zip.writestr("Metabolite (Organ) only.csv", 
                            Num[(Num['OrganTran_number'] ==0) & (Num['BloodMeta_number'] == 0)& (Num['OrganMeta_number'] > 0)].to_csv(index=False))                   
        with open("Transporter.zip", "rb") as file: 
            st.download_button(label = "Download transporter data",data = file,file_name = "Transporter.zip")
        del ALL
        gc.collect()
        #Upload
        degree1 = list(dict(nx.degree(G1)).values())
        x1 = [i for i in range(max(degree1)+1)]
        y1 = [degree1.count(i) for i in x1]
        Pk10=pd.DataFrame([x1,y1]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk11=Pk10[Pk10["k"]>0]
        Pk12=Pk11[Pk11["P(k)"]>0]

        #Database
        CPD1=pd.read_csv("./Database/230228CPD2Transporter.csv")
        CPD1["CPD"]='Organ_' + CPD1['CPD'].astype(str)
        CPD1=pd.concat([pd.read_csv("./Database/230228CPD2Transporter.csv"),CPD1])
        CPD1=CPD1.rename(columns={'CPD': 'Molecule'})
        CPD1=pd.concat([CPD1,ENSMUSG])

        Gd = nx.Graph()
        Gd.add_nodes_from(CPD1["Molecule"].unique())
        Gd.add_nodes_from(CPD1["Transporter"].unique())
        edge_lists=[tuple(x) for x in CPD1[["Molecule","Transporter"]].values]
        Gd.add_edges_from(edge_lists)
        degree = list(dict(nx.degree(Gd)).values())
        x = [i for i in range(max(degree)+1)]
        y = [degree.count(i) for i in x]
        Pk=pd.DataFrame([x,y]).T.rename(columns={0: 'k'}).rename(columns={1: 'P(k)'})
        Pk1=Pk[Pk["k"]>0]
        Pk2=Pk1[Pk1["P(k)"]>0]

        st.write(pd.DataFrame(data=np.array([[nx.density(Gd),nx.degree_pearson_correlation_coefficient(Gd,0),np.mean(degree)],
                        [nx.density(G1),nx.degree_pearson_correlation_coefficient(G1,0),np.mean(degree1)]]),
                    index=['Database', 'Uploaded data'],columns=['Density', 'Assortativity','Mean degree']))
        f=open('./Fig/D4.txt', 'r')
        st.write(f.read())   
        
        fig, ax = plt.subplots(figsize=(2, 2))
        #Dat
        df_X =np.log10(Pk2["k"])
        df_y =np.log10(Pk2["P(k)"])
        color="gray"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Database:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")

        #Upload
        df_X =np.log10(Pk12["k"])
        df_y =np.log10(Pk12["P(k)"])
        color="black"
        ax.scatter(df_X, df_y,c=color,s=1)
        x=np.arange(0, np.max(df_X),0.1)
        df_X = sm.add_constant(df_X)
        model = sm.OLS(df_y, df_X)
        result = model.fit()
        ax.plot(x,result.params["k"]*x+result.params["const"],c=color,
                label="Upload:"+str(-result.params["k"])[0:4]+" ("+str(-result.conf_int(alpha=0.05).T["k"][1])[0:4]+" - "+str(-result.conf_int(alpha=0.05).T["k"][0])[0:4]+")")



        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel("log"+"${_1}$"+"${_0}$"+"k")
        ax.set_ylabel("log"+"${_1}$"+"${_0}$"+"N(k)")
        ax.legend(bbox_to_anchor=(2.5, 0.5), loc="upper right",fontsize=8)
        st.pyplot(fig) 
        
        f=open('./Fig/D5.txt', 'r')
        st.write(f.read())
        del CPD1
        del CPDA1
        gc.collect()
        #st.stop()
        
        current_variables = list(globals().keys())
        exclude_list = ['current_variables', 'exclude_list','selected_option']
        variables_to_delete = [var for var in current_variables if var not in exclude_list]

        for var_name in variables_to_delete:
            del globals()[var_name]
        import gc
        gc.collect()                         
    else:
        st.write("Please upload metabolome (organ or cell), metabolome (blood or medium) and transciptome data (organ or cell).")
        