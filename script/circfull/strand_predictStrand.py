import numpy as np,pandas as pd,gzip,os,sys
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import roc_curve, auc, classification_report
from sklearn.ensemble import RandomForestClassifier

def predictStrand(outPrefix,fastqFile):
    consFL_strand=pd.read_csv(outPrefix+'consFL_strand.txt',sep='\t',low_memory=False)
    consFL_strand=consFL_strand.set_index('ID')
    consFL_strand.loc[consFL_strand.index,'rawID']=consFL_strand.index
    consFL_strand_counts=consFL_strand['rawID'].value_counts()
    consFL_strand_counts_uq=list(consFL_strand_counts[consFL_strand_counts>1].index)
    consFL_strand_PASS=consFL_strand.loc[[i not in consFL_strand_counts_uq for i in consFL_strand['rawID']]]
    consFL_strand_PASS=consFL_strand_PASS.dropna()
    consFL_strand_PASS=consFL_strand_PASS.loc[consFL_strand_PASS['exonNum']>1]
    strandSeq=consFL_strand_PASS.loc[:,['rawID','seqStrand']]
    strandSeq=strandSeq.set_index('rawID')
    strandDict={'+':1,'-':0}
    strandSeq['value']=strandSeq.apply(lambda x: strandDict[x['seqStrand']],axis=1)
    primerScore=pd.read_csv(outPrefix+'primer_score.txt',sep='\t',names=["readName","AnchorX","AnchorY","polyT","AnchorX_rev","AnchorY_rev","polyA"])
    primerScore.index=primerScore['readName']
    if os.path.exists(outPrefix+'rawReadLen.txt'):
        seqLenDF=pd.read_csv(outPrefix+'rawReadLen.txt',header=None,sep='\t',index_col=0)
    else:
        fastqFile=open(fastqFile)
        seqLen={}
        while True:
            fq1=fastqFile.readline()
            if fq1:
                fq2=fastqFile.readline().strip()
                fq3=fastqFile.readline()
                fq4=fastqFile.readline()
                ID=fq1.split()[0][1:]
                seqLen[ID]=len(fq2)
            else:
                break
        fastqFile.close()
        seqLenDF=pd.Series(seqLen)
        seqLenDF.to_csv(outPrefix+'rawReadLen.txt',sep='\t',header=True)
        
    primerScore['len']=seqLenDF.loc[primerScore.index]
    primer_len=[23,24,24,23,24,24]

    primerScore_mat=primerScore.loc[:,["AnchorX","AnchorY","polyT","AnchorX_rev","AnchorY_rev","polyA"]].copy()


    primerScore_mat=primerScore_mat/primer_len
    primerScore_mat['XY']=primerScore_mat['AnchorX']-primerScore_mat['AnchorY']
    primerScore_mat['XY_rev']=primerScore_mat['AnchorX_rev']-primerScore_mat['AnchorY_rev']
    primerScore_mat['XX']=primerScore_mat['AnchorX']-primerScore_mat['AnchorX_rev']
    primerScore_mat['YY']=primerScore_mat['AnchorY']-primerScore_mat['AnchorY_rev']
    primerScore_mat['TA']=primerScore_mat['polyT']-primerScore_mat['polyA']
    primerScore_mat['len']=primerScore['len']
    primerScore_mat=primerScore_mat.loc[[not i for i in  primerScore_mat.index.duplicated()]]

    primerScore_sub=primerScore_mat.loc[primerScore_mat.index.intersection(strandSeq.index),].reindex(strandSeq.index)
    primerScore_sub['value']=list(strandSeq['value'])
    primerScore_sub=primerScore_sub.dropna()

    X_train,X_test,Y_train,Y_test=train_test_split(primerScore_sub.loc[:,primerScore_sub.columns!='value'],primerScore_sub['value'],stratify=primerScore_sub['value'],random_state=66)

    training_accuracy = []
    test_accuracy = []
    rf=RandomForestClassifier(random_state=0,n_estimators=200,min_samples_leaf=1, min_samples_split=10)
    rf.fit(X_train,Y_train)
    print("Feature importance: ")
    print(rf.feature_importances_)
    print("Accuracy on training set: {:.3f}".format(rf.score(X_train,Y_train)))
    print("Accuracy on test set: {:.3f}".format(rf.score(X_test,Y_test)))
    prob_predict_y_validation = rf.predict_proba(X_test)
    predictions_validation = prob_predict_y_validation[:, 1]  
    fpr, tpr, _ = roc_curve(Y_test, predictions_validation)  
    roc_auc = auc(fpr, tpr)
    print("AUC: {:.3f}".format(roc_auc))
    '''
    Feature importance: 
    [0.16590163 0.02183294 0.0891461  0.02298425 0.01096239 0.15567186
     0.11377447 0.11286144 0.03223863 0.01393638 0.23858744 0.02210247]
    Accuracy on training set: 0.975
    Accuracy on test set: 0.951
    AUC: 0.988
    '''
    predictStrand=rf.predict_proba(primerScore_mat)
    predictStrand_sense=predictStrand[:,1] # the probability of reads have the same direction as RNA 
    predictStrandDf=pd.Series(predictStrand_sense,index=primerScore_mat.index)
    predictStrandDf.to_csv(outPrefix+'strandProbability.txt',sep='\t',header=False)