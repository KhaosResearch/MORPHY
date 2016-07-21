#include "Tree.h"

Tree::Tree() {
}

//Constructor that creates the created Tree based on the input Tree file and the Given Sequences
Tree::Tree(string filename, Sequences* Seqs) {
    string line;
    ifstream myfile (filename.c_str());

    if (myfile.is_open()) {
        while ( myfile.good() ) { //For each Node
            getline (myfile,line);
            if (line.compare("")==0 || line.compare("\r")==0  || line.compare("\r\n")==0 || line.compare("\n")==0 || line.compare(" ")==0) { continue; cout<<"Skipped"<<endl; } //skip empty lines

            string NodeLbl = "";
            double NodeLen = 0.0;

            if (line.find(":")!= string::npos) {
                int delimPos = line.find(":");
                NodeLbl = line.substr(0, delimPos);
                NodeLen = atof(line.substr(delimPos+1, line.size()-delimPos).c_str());
            }

            if ( GetPosOfElement(UniqueLengths, NodeLen, 0) == -1 ) {
                UniqueLengths.push_back(NodeLen);
            }


            if (NodeLbl.compare("A")==0 || NodeLbl.compare("B")==0 || NodeLbl.compare("C")==0 || NodeLbl.compare("D")==0){
                Node temp;
                temp.N1.Length = NodeLen;
                temp.N2.Length = -1.0f;
                temp.N3.Length = -1.0f;
                TreeNodes.push_back(temp);
            }else if (NodeLbl.compare("E")==0){
                Node temp;
                temp.N1.Length = NodeLen;
                temp.N2.Length = TreeNodes[0].N1.Length;
                temp.N3.Length = TreeNodes[1].N1.Length;
                TreeNodes.push_back(temp);
            }else if (NodeLbl.compare("F")==0){
                Node temp;
                temp.N1.Length = NodeLen;
                temp.N2.Length = TreeNodes[2].N1.Length;
                temp.N3.Length = TreeNodes[3].N1.Length;
                TreeNodes.push_back(temp);
            }
        }
    } else cerr << "Unable to open file: "<<filename<<" for reading.";
    myfile.close();
}

Tree::~Tree() {
    TreeNodes.clear();
    UniqueLengths.clear();
}

double Tree::GetLengthOfNode(int NodeIdx, int InnerNodeIdx) {
    switch(InnerNodeIdx){
    case 1: return TreeNodes[NodeIdx].N1.Length;
    case 2: return TreeNodes[NodeIdx].N2.Length;
    case 3: return TreeNodes[NodeIdx].N3.Length;
    }
}

DoubleMat Tree::GetFPTableOfNode(int NodeIdx, int InnerNodeIdx) {
    switch(InnerNodeIdx){
    case 1: return TreeNodes[NodeIdx].N1.FPTable;
    case 2: return TreeNodes[NodeIdx].N2.FPTable;
    case 3: return TreeNodes[NodeIdx].N3.FPTable;
    }

}

//Definition of the toString Operation
ostream& operator<<(ostream& ostr, const Tree& tree) {
    for (int i=0;i<tree.Get_TreeNodes().size();i++){
        if      (i==0) { ostr<<"Node A:\n"; }
        else if (i==1) { ostr<<"Node B:\n"; }
        else if (i==2) { ostr<<"Node C:\n"; }
        else if (i==3) { ostr<<"Node D:\n"; }
        else if (i==4) { ostr<<"Node F:\n"; }
        else if (i==5) { ostr<<"Node E:\n"; }
        ostr<<"\tInnerNode N1 (Length:" <<tree.Get_TreeNodes()[i].N1.Length << ")"<<endl;
        for (int k=0;k<tree.Get_TreeNodes()[i].N1.FPTable.size();k++){
            ostr<<"\t  ";
            for (int l=0;l<tree.Get_TreeNodes()[i].N1.FPTable[k].size();l++){
                ostr<<setw(10)<<tree.Get_TreeNodes()[i].N1.FPTable[k][l];
            }
            ostr<<endl;
        }
        ostr<<endl;

        ostr<<"\tInnerNode N2 (Length:" <<tree.Get_TreeNodes()[i].N2.Length << ")"<<endl;
        for (int k=0;k<tree.Get_TreeNodes()[i].N2.FPTable.size();k++){
            ostr<<"\t  ";
            for (int l=0;l<tree.Get_TreeNodes()[i].N2.FPTable[k].size();l++){
                ostr<<setw(10)<<tree.Get_TreeNodes()[i].N2.FPTable[k][l];
            }
            ostr<<endl;
        }
        ostr<<endl;


        ostr<<"\tInnerNode N3 (Length:" <<tree.Get_TreeNodes()[i].N3.Length << ")"<<endl;
        for (int k=0;k<tree.Get_TreeNodes()[i].N3.FPTable.size();k++){
            ostr<<"\t  ";
            for (int l=0;l<tree.Get_TreeNodes()[i].N3.FPTable[k].size();l++){
                ostr<<setw(10)<<tree.Get_TreeNodes()[i].N3.FPTable[k][l];
            }
            ostr<<endl;
        }
        ostr<<endl;

    }
    return ostr;
}
