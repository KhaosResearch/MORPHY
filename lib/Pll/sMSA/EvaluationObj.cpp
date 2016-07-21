#include "EvaluationObj.h"

//Default Constructor
EvaluationObj::EvaluationObj()  {
}

//Virtual Deconstructor
EvaluationObj::~EvaluationObj() { }

//Function to convert a Given alignment into a matrix of indices representing the alignment
void EvaluationObj::ConvertAlignmentToMatrix(StringVecRef alignment, IntMatRef mat) {
    mat.resize(alignment.size(), IntVec(alignment[0].size(), -1));

    for (int i=0; i<alignment.size(); i++) {
        int counter = 0;
        for (int j=0; j<alignment[i].size(); j++) {
            if (alignment[i][j]!='-') {
                mat[i][j] = counter++;
            }
        }
    }
}

//Function to return the MCC score after comparing the predicted and the reference alignment at a pairwise level
string EvaluationObj::GetMCCScore(IntMatRef predMat, IntMatRef refMat) {
    int TP = 0, FP = 0, TN = 0, FN = 0;
    //For every position in the alignment matrix..
    //omp_set_num_threads( 4 );
    #pragma omp parallel for
    for (int k=0; k<predMat[0].size(); k++) {
        //Compare all pairs between predicted and reference matrices
        for (int i=0; i<predMat.size()-1; i++) { //last sequence has no point being compared
            for (int j=(i+1); j<predMat.size(); j++) {
                if (predMat[i][k]==-1 && predMat[j][k]==-1) continue;  //Ignore gapped positions in the prediction column

                int l, m;
                if (predMat[i][k]!=-1) l = GetPosOfElement(refMat[i], predMat[i][k], 0);
                if (predMat[j][k]!=-1) m = GetPosOfElement(refMat[j], predMat[j][k], 0);

                if (refMat[i][l]==-1 && refMat[j][m]==-1) continue;
/*
                if (predMat[j][k]!=-1 && refMat[j][p]!=-1 && predMat[j][k]==refMat[j][p]) {
                    TP++;
                } else if (predMat[j][k]!=-1 && refMat[j][p]==-1) {
                    FP++;
                } else if (predMat[j][k]==-1 && refMat[j][p]==-1) {
                    TN++;
                } else if (predMat[j][k]==-1 && refMat[j][p]!=-1) {
                    FN++;
                } else { //i.e. predMat[j][k]!=-1 && refMat[j][p]!=-1 && predMat[j][k]!=refMat[j][p]
                    FP++;
                }
*/
            }
        }
    };

    MCCScoreSummary score;
    if ((TP+FN)!=0) {
        score.SensitivityScore = (double) TP/(double)(TP+FN);
    } else{
        score.SensitivityScore = (double) TP;
    }

    if ((TP+FP)!=0) {
        score.SelectivityScore = (double) TP/(double)(TP+FP);
    } else{
        score.SelectivityScore = (double) TP;
    }

    double nominator   = (double) pow(TP,2);
    double denominator = (double)((TP+FN)*(TP+FP));
    if (denominator==0) {
        score.MCCScoreValue = sqrt(nominator);//artitrarily set denominator equal to 1
    } else {
        score.MCCScoreValue = sqrt(nominator/denominator);
    }

    stringstream ss; ss.str("");
    ss<<"                      MCC score: "<<score.MCCScoreValue<<endl;
    ss<<"         Sensitivity score (TP/(TP+FN)): "<<score.SensitivityScore<<endl;
    ss<<"         Selectivity score (TP/(TP+FP)): "<<score.SelectivityScore<<endl;
    ss<<endl;
    return ss.str();
}

//Function to return the similarity-based SP (Sum of Pairs) score at a pairwise level. Shows how close the aligned sequences are to each other
string EvaluationObj::GetSPSimScore(StringVecRef align){
    int match = 1   , mismatch   = 0;
    int simScore = 0, totalPairs = 0;

    //For every position in the alignment matrix..
    //omp_set_num_threads( 4 );
    //#pragma omp parallel for
    for (int k=0; k<align[0].size(); k++) {
        for (int i=0;i<align.size()-1;i++){
            //Neglect Gaps
            if (align[i][k]=='-'){
                continue;
            }
            for (int j=(i+1); j<align.size(); j++){
                //Neglect Gaps
                if (align[j][k]=='-'){
                    continue;
                }

                if (align[i][k]==align[j][k]) {
                    simScore += match;
                }else{
                    simScore += mismatch;
                }
                totalPairs++;
            }
        }
    }

    SPSimScoreSummary score;
    score.SPSimScore = simScore;
    score.TotalNumOfPairs = totalPairs;
    score.SPScoreValue = (double)simScore/(double)totalPairs;


    stringstream ss; ss.str("");
    ss<<"               Alignment (Sim) SP Score: "<<left<<setw(9)<<score.SPScoreValue<<left<<" ("<<score.SPSimScore<<"/"<<score.TotalNumOfPairs<<")"<<endl;

    return ss.str();
}

//Function to return the distance-based SP (Sum of Pairs) score at a pairwise level. Shows how distant the aligned sequences are to each other
string EvaluationObj::GetSPDistScore(StringVecRef align){
    int match = 0,     mismatch = 1;
    int distScore = 0, totalPairs = 0;

    //For every position in the alignment matrix..
    //omp_set_num_threads( 4 );
    //#pragma omp parallel for
    for (int k=0; k<align[0].size(); k++) {
        for (int i=0;i<align.size()-1;i++){
            //Neglect Gaps
            if (align[i][k]=='-'){
                continue;
            }
            for (int j=(i+1); j<align.size(); j++){
                //Neglect Gaps
                if (align[j][k]=='-'){
                    continue;
                }

                if (align[i][k]==align[j][k]) {
                    distScore += match;
                }else{
                    distScore += mismatch;
                }
                totalPairs++;
            }
        }
    }

    SPDistScoreSummary score;
    score.SPDistScore = distScore;
    score.TotalNumOfPairs = totalPairs;
    score.SPScoreValue = (double)distScore/(double)totalPairs;

    stringstream ss; ss.str("");
    ss<<"               Alignment(Dist) SP Score: "<<left<<setw(9)<<score.SPScoreValue<<left<<" ("<<score.SPDistScore<<"/"<<score.TotalNumOfPairs<<")"<<endl;

    return ss.str();

}

//Function to return the TC score after comparing the predicted and the reference alignment at a columnwise level
string EvaluationObj::GetTCScore(StringVecRef pred, StringVecRef ref) {
    TCScoreSummary score;
    score.RefColumns = ref[0].size();

    int numOfMatches = 0;

    IntHyp Idxs(pred.size(), IntMat(5, IntVec(0,0)));
    //#pragma omp parallel for
    for (int i=0;i<pred.size();i++){
        for (int j=0;j<pred[i].size();j++){
            if       (pred[i][j]=='A'){
                Idxs[i][0].push_back(j);
            }else if (pred[i][j]=='C'){
                Idxs[i][1].push_back(j);
            }else if (pred[i][j]=='G'){
                Idxs[i][2].push_back(j);
            }else if (pred[i][j]=='T'){
                Idxs[i][3].push_back(j);
            }else if (pred[i][j]=='-'){
                Idxs[i][4].push_back(j);
            }
        }
    }

    //#pragma omp parallel for
    for (int j=0;j<ref[0].size();j++){
        vector< vector<int>* > CompMat;
        for (int i=0;i<ref.size();i++) {
            IntVecPtr Tail;
            if      (ref[i][j]=='A') { Tail = &(Idxs[i][0]); }
            else if (ref[i][j]=='C') { Tail = &(Idxs[i][1]); }
            else if (ref[i][j]=='G') { Tail = &(Idxs[i][2]); }
            else if (ref[i][j]=='T') { Tail = &(Idxs[i][3]); }
            else if (ref[i][j]=='-') { Tail = &(Idxs[i][4]); }
            CompMat.push_back(Tail);
        }


        for (int k=0;k<(*(CompMat[0])).size();k++) {
            int counter = 1;
            IntVec DelPos;
            DelPos.push_back(k);
            bool matchFound = false;

            for (int i=1;i<CompMat.size();i++) {
                int p = GetPosOfElement((*(CompMat[i])), (*(CompMat[0]))[k], 0);
                if (p==-1) { break;     }
                else       { counter++; DelPos.push_back(p); }
            }

            if (counter==CompMat.size()){
                matchFound = true;
                numOfMatches++;

                for (int l=0;l<DelPos.size();l++){
                    (*(CompMat[l])).erase((*(CompMat[l])).begin() + DelPos[l]);
                }

                break;
            }
        }
    }

    score.SharedColumns = numOfMatches;
    score.TCScoreValue = (double) numOfMatches/ (double) ref[0].size();

    stringstream ss; ss.str("");
    ss<<"                   Column-wise TC Score: "<<left<<setw(9)<<score.TCScoreValue<<left<<" ("<<score.SharedColumns<<"/"<<score.RefColumns<<")"<<endl;

    return ss.str();
}
