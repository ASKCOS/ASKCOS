import numpy as np
from sklearn.metrics.ranking import auc

def get_AUROC(true_val, pred_val):
    '''
    data should contain following:
    list of true classifications 
    list of predicted values
    '''
    
    TPR = [1.0]
    FPR = [1.0]
    
    for t in range(1,999):
        TP = 0.0
        FP = 0.0
        FN = 0.0
        TN = 0.0
        threshold = float(t)/1000.0
        for i,data_point in enumerate(pred_val):
            if data_point > threshold:
                if true_val[i] == 1:
                    TP += 1.0
                else:
                    FP += 1.0
            else:
                if true_val[i] == 1:
                    FN += 1.0
                else:
                    TN += 1.0
        if TP+FN == 0 or FP+TN == 0:
            continue
        TPR_i = TP/(TP+FN)
        FPR_i = FP/(FP+TN)
        
        TPR.append(TPR_i)
        FPR.append(FPR_i)
    TPR.append(0.0)
    FPR.append(0.0)
    auc = np.trapz(TPR, FPR)
    return auc, TPR, FPR

if __name__ == '__main__':
    true = [1,0,0,0,1,1,1,0,1,0,0,1]
    pred = [0.9, 0.1,0.2,0.6,0.7,0.8,0.5,0.7,0.4,0.3,0.9]
    (auc, TPR, FPR) = get_AUROC(true, pred)
    print auc
    print TPR
    print FPR
                    
            