#include<iostream>
#include<vector>
#include<cmath>
#include<map>
#include<time.h>
#include<fstream>
#include<algorithm>

using namespace std;

#include<iostream>
#include<vector>
#include<cmath>
#include<cstdlib>

#define UINT32 unsigned int
#define REAL double

//#define DEBUG

#define GP_AND 0
#define GP_ADD 1
#define GP_OR  2
#define GP_SUB 3
#define GP_XOR 4
#define GP_NOT 5
#define GP_SR1 6
#define GP_SL1 7

using namespace std;

typedef void (*gp_function)(vector<UINT32>&, vector<UINT32>&, vector<UINT32>&);

void gp_add(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
        r[i] = (a[i] + b[i]) & 15;
    }
};

void gp_sub(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
		r[i] = (a[i] - b[i]) & 15;
    }
};
void gp_and(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
		r[i] = (a[i] & b[i]) & 15;
    }
};
void gp_or(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
		r[i] = (a[i] | b[i]) & 15;
    }
};

void gp_xor(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
		r[i] = (a[i] ^ b[i]) & 15;
    }
};

void gp_not(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
        r[i] = ~a[i] & 15;
    }
};

void gp_sr1(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
        r[i] = (a[i] >> 1) & 15;
    }
};

void gp_sl1(vector<UINT32> &r, vector<UINT32> &a, vector<UINT32> &b){
    for(int i=0;i<(int)r.size();i++){
        r[i] = (a[i] << 1) & 15;
    }
};

int gp_evaluation(vector<UINT32> &r, vector<UINT32> &p){
    int sum = 0;
	for(int i=0;i<(int)r.size();i++){
        sum += abs((int)r[i] - (int)p[i]);
    }
    return sum;
};


class Symbol{
private:
    string symbol;
    int    numChildren;
public:
    gp_function    func;
    vector<UINT32> value;
    Symbol(string in_symbol, gp_function gf, int nc){
        symbol = in_symbol;
        func   = gf;
        numChildren = nc;
    };
    Symbol(string line){
        vector<string> s_gt;
        char *endp;
        numChildren = 0;
        func = NULL;
        s_gt = split(line, ",");
        symbol = s_gt[0];
        value.resize((int)(s_gt.size()-1));
        for(int i=1;i<(int)s_gt.size();i++){
            value[i-1] = strtoul(s_gt[i].c_str(), &endp, 10);
        }
    };
    ~Symbol(){
    };
    vector<string> split(string str, string delim) {
        vector<string> items;
        size_t dlm_idx;
        if(str.npos == (dlm_idx = str.find_first_of(delim))) {
            items.push_back(str.substr(0, dlm_idx));
        }
        while(str.npos != (dlm_idx = str.find_first_of(delim))) {
            if(str.npos == str.find_first_not_of(delim)) {
                break;
            }
            items.push_back(str.substr(0, dlm_idx));
            dlm_idx++;
            str = str.erase(0, dlm_idx);
            if(str.npos == str.find_first_of(delim) && "" != str) {
                items.push_back(str);
                break;
            }
        }
        return items;
    }
    void setValue(vector<UINT32> in_value){
        value.assign((int)in_value.size(), 0);
        for(int i=0;i<(int)value.size();i++){
            value[i] = in_value[i];
        }
    }
    void setSymbol(string a){symbol = a;};
    void setFunc(gp_function in_func){func = in_func;};
    gp_function getFunc(){return func;};
    string getSymbol(){return symbol;};
    UINT32 getValue(int idx){return value[idx];};
    int getValueLength(){return value.size();};
    int getNumChildren(){return numChildren;};
};

class BinaryTree{
private:
    map<int, Symbol*> tree;
    int maxDepth;
    vector<Symbol*> gpFunctions;
    vector<Symbol*> gpTerminals;
    Symbol* gpAnswer;
    UINT32 score;
public:
    BinaryTree(int maxdepth_in,
               vector<Symbol*>& gp_funcs,
               vector<Symbol*>& gp_terms,
               Symbol* gp_answer){
        maxDepth = maxdepth_in;
        gpFunctions = gp_funcs;
        gpTerminals = gp_terms;
        gpAnswer = gp_answer;
    };
    static bool cmpTree(BinaryTree* lhs, BinaryTree* rhs){
        bool t;
        // スコアが小さいほうが先
        if(lhs->score < rhs->score){
            t = true;
        }
        else{
            // スコアが等しいとき
            if(lhs->score == rhs->score){
                if(lhs->score == 0 && rhs->score==0){
                    if(lhs->countNode() < rhs->countNode()){
                        t = true;
                    }
                    else{
                        t = false;
                    }
                }
                else{
                    if(lhs->countNode() > rhs->countNode()){
                        t = true;
                    }
                    else{
                        if(lhs->countNode() == rhs->countNode()){
                            if(lhs->traverse() > rhs->traverse()){
                                t = true;
                            }
                            else{
                                t = false;
                            }
                        }
                        else{
                            t = false;
                        }
                    }
                }
            }
            // スコアが大きい
            else{
                t = false;
            }
        }
        return t;
    };

    static bool cmpTree2(BinaryTree* lhs, BinaryTree* rhs){
        bool t;
        // スコアが小さいほうが先
        if(lhs->score < rhs->score){
            t = true;
        }
        else{
            // スコアが等しいとき
            if(lhs->score == rhs->score){
                if(lhs->score == 0 && rhs->score==0){
                    if(lhs->countNode() < rhs->countNode()){
                        t = true;
                    }
                    else{
                        t = false;
                    }
                }
                else{
                    if(lhs->countNode() < rhs->countNode()){
                        t = true;
                    }
                    else{
                        if(lhs->countNode() == rhs->countNode()){
                            if(lhs->traverse() > rhs->traverse()){
                                t = true;
                            }
                            else{
                                t = false;
                            }
                        }
                        else{
                            t = false;
                        }
                    }
                }
            }
            // スコアが大きい
            else{
                t = false;
            }
        }
        return t;
    };

    int countNode(){ return (int)tree.size(); };
    void generateTree(int idx);
    map<int, Symbol*>& getTree(){return tree;};
    int  selectNodeKey(int n);
    int getDepth(int idx);
    void computeScore();
    void setScore(int score);
    UINT32 getScore(){ return score; };
    int  getNodeNum(){return (int)tree.size();};
    Symbol* getAnswer(){return gpAnswer;};
    string traverse();
    void mutation();
    int  getHeight(int idx);
    void deleteNode(int idx);
    BinaryTree* clone();
    vector<UINT32> calculateScore(int idx);
};

/*
    map<int, Symbol*> tree;
    int maxDepth;
    vector<Symbol*> gpFunctions;
    vector<Symbol*> gpTerminals;
    Symbol* gpAnswer;
    UINT32 score;
*/

int BinaryTree::getDepth(int idx){
    int n;
    n = 0;
    while((int)pow(2,n) < idx){
        n++;
    }
    return n;
};


void BinaryTree::computeScore(){
    vector<UINT32> tree_output;
    tree_output = calculateScore(1);
    score = gp_evaluation(gpAnswer->value, tree_output);
};

void BinaryTree::setScore(int in_score){
    score = in_score;
};

BinaryTree* BinaryTree::clone(){
    BinaryTree* bt;
    bt = new BinaryTree(maxDepth, gpFunctions, gpTerminals, gpAnswer);
    // getTree
    bt->tree = tree;
    bt->setScore(score);
    /*
    cout << "*==== clone ====*" << endl;
    cout << tree.size() << endl;
    cout << bt->tree.size() << endl;
    */
    return bt;
};

// idx以下のノードの高さを計算する(1オリジン)
int BinaryTree::getHeight(int idx){
    int num = 0;
    if(tree.find(idx) != tree.end()){
        if(getHeight(2*idx+0) < getHeight(2*idx+1))
            num = getHeight(2*idx+1) + 1;
        else
            num = getHeight(2*idx+0) + 1;
    }
    else{
        num = 0;
    }
    return num;
};

vector<UINT32> BinaryTree::calculateScore(int idx){
    Symbol* s;
    vector<UINT32> ret;
    vector<UINT32> u1, u2;
    // もしも終端ノードであれば数値を取り出す

    // cout << "Calculate Score :[" << idx << "]" << endl;
    // キーidxの存在を確認
    if(tree.find(idx) != tree.end()){
        s = tree[idx];
        // cout << s->getSymbol() << endl;
        if(s->getNumChildren() == 0){
            ret.resize(s->getValueLength());
            for(int i=0;i<s->getValueLength();i++){
                ret[i] = s->getValue(i);
            }
        }
        // 終端ノードでなければ計算する
        else{
            u1 = calculateScore(2*idx+0);
            u2 = calculateScore(2*idx+1);
            ret.resize((int)u1.size());
            s->func(ret, u1, u2);
        }
    }
    return ret;
};

// idx以下の木を作る．以下の木は存在しないことが前提
void BinaryTree::generateTree(int idx){
    int r;
    Symbol* s;

    // cout << "generateNode [" << idx << "]" << endl;
    r = rand() % RAND_MAX;
    // 子ノードが作れる位置の場合は，終端または演算ノードを置く
    if(idx < pow(2, maxDepth-1)){
        r = r % ((int)gpTerminals.size() + (int)gpFunctions.size());
        if(r < (int)gpTerminals.size()){
            tree.insert(make_pair(idx, gpTerminals[r]));
        }
        else{
            r -= (int)gpTerminals.size();
            tree.insert(make_pair(idx, gpFunctions[r]));
            for(int i=0;i<gpFunctions[r]->getNumChildren();i++){
                generateTree(2*idx+i);
            }
        }
    }
    else{
        // 端に到達している場合は，終端ノードだけ作る
        r = r % (int)gpTerminals.size();
        tree.insert(make_pair(idx, gpTerminals[r]));
    }
};

// idx以下の木を削除する
void BinaryTree::deleteNode(int idx){
    map<int, Symbol*>::iterator it;
    it = tree.find(idx);
    if(it != tree.end()){
        tree.erase(idx);
        // cout << "delete t["<< idx << "]" << endl;
        deleteNode(2*idx+0);
        deleteNode(2*idx+1);
    }
};

// 幅優先探索
string BinaryTree::traverse(){
    string s;
    map<int, Symbol*>::iterator it;
    it = tree.begin();
    /*
    while(it != tree.end()){
        // cout << it->first << " : " << it->second->getSymbol() << endl;
        cout << "{" << it->first << ", " << it->second->getSymbol() << "} ";
        it++;
    }
    */
    while(it != tree.end()){
        s += it->second->getSymbol() + ",";
        it++;
    }
    return s;
};

int BinaryTree::selectNodeKey(int n){
    int idx;
    map<int, Symbol*>::iterator it;
    it = tree.begin();
    advance(it, n);
    // cout << n << " : " << it->first << endl;
    idx = it->first;
    return idx;
};


void BinaryTree::mutation(){
    int ir, idx;
    REAL rr;
    map<int, Symbol*>::iterator it;
    
    // rrで突然変異発生確率
    rr = rand();
    // r = (REAL)rand()/((REAL)RAND_MAX+1.0) * (REAL)sum;

    // [1, tree.size()]までの乱数ir生成
	// tree.at(ir)でir番目の木(のvalue)が取り出せる
    ir = (rand() % RAND_MAX) % tree.size();
    idx = selectNodeKey(ir);
    /*
    cout << " ==== Mutation ==== " << endl;
    cout << "  node size:" << tree.size() << endl;
    cout << "  [" << ir << "] = t[" << idx << "] is mutated." << endl;
    cout << "  Before: " << endl;
    traverse();
    */
    deleteNode(idx);
    generateTree(idx);
    /*
    cout << "  After: " << endl;
    traverse();
    */
    // cout << " ==== End Mutation ==== " << endl;
};

class Population{
private:
    vector<BinaryTree*> individual;
    vector<BinaryTree*> elite;
    int maxDepth;
public:
    BinaryTree* getIndividual(int idx){ return individual[idx]; };
    int getNumIndividual(){ return (int)individual.size(); };
    int getNumElite(){ return (int)elite.size(); };
    BinaryTree* getElite(int idx){ return elite[idx]; };
    Population(int btnum, int elnum, int maxDepth,
               vector<Symbol*>& gpFunctions,
               vector<Symbol*>& gpTerminals,
               Symbol* gpAnswer);
    void preserveElite(map<string, BinaryTree*>& repo, int flag);
    // void preserveElite();
    void selection();
    void mutation();
    void i_pop_back(){individual.pop_back();};
    void i_push_back(BinaryTree* t){individual.push_back(t);};
    void crossOver();
	void swapNode(map<int, Symbol*>& bt1, map<int, Symbol*>& bt2,
                  int idx1, int idx2);
};

Population::Population(int btnum, int elnum, int in_maxDepth,
                       vector<Symbol*>& gpFunctions,
                       vector<Symbol*>& gpTerminals,
                       Symbol* gpAnswer){
    //individual = new BinaryTree(maxDepth, gpFunctions, gpTerminals);
    maxDepth = in_maxDepth;
    individual.resize(btnum);
    cout << " btnum:" << btnum << endl;
    for(int i=0;i<btnum;i++){
        individual[i] = new BinaryTree(maxDepth, gpFunctions,
                                       gpTerminals, gpAnswer);
        individual[i]->generateTree(1);
    }
    for(int i=0;i<(int)individual.size();i++){
        cout << individual[i]->traverse() << endl;
        individual[i]->computeScore();
    }
    elite.resize(elnum);
    for(int i=0;i<elnum;i++){
        elite[i] = individual[i]->clone();
    }
}

void Population::preserveElite(map<string, BinaryTree*>& repo, int flag){
    vector<BinaryTree*> p;
    map<string, BinaryTree*> m;
    int elnum, pnum, j,k;
    // cout << " == S preserveElite ==" << endl;
    p.clear();
    p.resize((int)individual.size()+(int)elite.size());
    for(int i=0;i<(int)elite.size();i++){
        p[i] = elite[i];
    }
    for(int i=0;i<(int)individual.size();i++){
        p[i+(int)elite.size()] = individual[i];
    }
    for(int i=0;i<(int)p.size();i++){
        p[i]->computeScore();
    }
    /*
    cout << "  == cs preserveElite ==" << endl;
    cout << "psize : " << p.size() << endl;
    */
    if(flag == 0){
        sort(p.begin(), p.end(), BinaryTree::cmpTree);
    }
    else{
        sort(p.begin(), p.end(), BinaryTree::cmpTree2);
    }
    //sort(individual.begin(), individual.end(), BinaryTree::cmpTree);
    elnum = (int)elite.size();
    m.clear();
    j=0;
    for(int i=0;i<(int)individual.size();i++){
        if(j < elnum){
            if((m.find(p[i]->traverse()) == m.end())
                && repo.find(p[i]->traverse()) == repo.end()){
                elite[j] = p[i]->clone();
                m.insert(make_pair(p[i]->traverse(), p[i]));
                repo.insert(make_pair(p[i]->traverse(), p[i]));
                j++;
            }
        }
    }
    k=0;
    while(j < elnum){
        elite[j] = p[k]->clone();
        j++; k++;
    }
    // cout << "  == cloneelite preserveElite ==" << endl;
    pnum = (int)p.size();
    for(int i=0;i<pnum;i++){
        if(i<(int)individual.size()){
            individual[i] = p[i];
        }
        else{
            delete p[i];
        }
    }
    p.clear();
    #ifdef DEBUG
    cout << " == E preserveElite ==" << endl;
    #endif
};

void Population::selection(){
    vector<BinaryTree*> nt;
    int max_score, idx;
    vector<REAL> iscore;
    REAL sum, r, rs, num;

    #ifdef DEBUG
    cout << " ==== Begin Selection ==== " << endl;
    #endif
    /*
    for(int i=0;i<(int)individual.size();i++){ individual[i]->computeScore(); }
    sort(individual.begin(), individual.end(), BinaryTree::cmpTree);
    */

    max_score = individual[(int)individual.size()-1]->getScore();
    iscore.resize((int)individual.size());
    sum = 0.0;
    for(int i=0;i<(int)individual.size();i++){
        iscore[i] = (REAL)max_score / ((REAL)individual[i]->getScore()+1.0);
        sum += iscore[i];
    }

    #ifdef DEBUG
    for(int i=0;i<(int)individual.size();i++){
        cout << i << " : " << individual[i]->getScore()
             << " : " << iscore[i] << endl;
    }
    cout << sum << " : " << max_score << endl;
    #endif
    
    nt.resize((int)individual.size());
    for(int i=0;i<(int)individual.size();i++){
        r = (REAL)rand()/((REAL)RAND_MAX+1.0) * (REAL)sum;
        idx = 0; rs = 0.0;
        while(rs < r){
            rs += iscore[idx];
            idx++;
        }
        idx -= 1;
        nt[i] = individual[idx]->clone();
    }
    
    num = (int)individual.size();
    for(int i=0;i<num;i++){
        delete individual[i];
        individual[i]=nt[i];
    }
    nt.clear();

};

void Population::swapNode(map<int, Symbol*>& bt1, map<int, Symbol*>& bt2,
                          int idx1, int idx2){
    Symbol *s1, *s2;

    // 木bt1のノードidx1が存在する
    if(bt1.find(idx1) != bt1.end()){
        // 木bt2のノードidx2も存在する
        if(bt2.find(idx2) != bt2.end()){
            s1 = bt1[idx1];
            bt1.erase(idx1);
            s2 = bt2[idx2];
            bt2.erase(idx2);
            bt1.insert(make_pair(idx1, s2));
            bt2.insert(make_pair(idx2, s1));
            swapNode(bt1, bt2, 2*idx1+0, 2*idx2+0);
            swapNode(bt1, bt2, 2*idx1+1, 2*idx2+1);
        }
        else{
            s1 = bt1[idx1];
            bt1.erase(idx1);
            bt2.insert(make_pair(idx2, s1));
            swapNode(bt1, bt2, 2*idx1+0, 2*idx2+0);
            swapNode(bt1, bt2, 2*idx1+1, 2*idx2+1);
        }
    }
    else{
        // bt1にidx1がなくて，bt2にidx2がある場合
        if(bt2.find(idx2) != bt2.end()){
            s2 = bt2[idx2];
            bt2.erase(idx2);
            bt1.insert(make_pair(idx1, s2));
            swapNode(bt1, bt2, 2*idx1+0, 2*idx2+0);
            swapNode(bt1, bt2, 2*idx1+1, 2*idx2+1);
        }
        else{
            // 両方存在しない場合は何もしない
        }
    }
};

void Population::mutation(){
    for(int i=0;i<(int)individual.size();i++){
        individual[i]->mutation();
    }
};

void Population::crossOver(){
    BinaryTree *bt1, *bt2;
    int idx1, idx2;
    for(int i=0;i<(int)individual.size()/2;i++){
        bt1 = individual[2*i+0];
        bt2 = individual[2*i+1];
        // 交叉対象となるノードを選択
        // (木1におけるidx1の現在深さ+木2におけるidx2以下の高さ) < 最大深さ
        // (木2におけるidx2の現在深さ+木1におけるidx1以下の高さ) < 最大深さ
        do{
            idx1 = bt1->selectNodeKey((rand() % RAND_MAX) % bt1->getNodeNum());
            idx2 = bt2->selectNodeKey((rand() % RAND_MAX) % bt2->getNodeNum());
        }while((bt1->getDepth(idx1)+bt2->getHeight(idx2) > maxDepth) ||
               (bt2->getDepth(idx2)+bt1->getHeight(idx1) > maxDepth));
        /*
        cout << "  bt1[" << 2*i+0 << "] key:" << idx1 << ", "
             << "  bt2[" << 2*i+1 << "] key:" << idx2 << endl;
        */
        swapNode(bt1->getTree(), bt2->getTree(), idx1, idx2);
    }
};

int main(int argc, char** argv){
    ifstream ifs;
    string buf;
    vector<Symbol*> gpFunctions;
    vector<Symbol*> gpTerminals;
    Symbol* gpAnswer;

	int cycleLimit    = 5000;
    int numPopulation = 8;
    int numIndividual = 1280 / numPopulation;
    int numElite      = 2;
    int maxDepth      = 12;
    int numMigration  = 1;
    int dispInterval   = 5;
    int migrateInterval = 20;
    vector<Population*> p;
    map<string, BinaryTree*> m;

    srand((UINT32)time(NULL));

	gpFunctions.resize(8);
    gpFunctions[GP_AND] = new Symbol("and", gp_and, 2);
    gpFunctions[GP_ADD] = new Symbol("add", gp_add, 2);
    gpFunctions[GP_OR]  = new Symbol("or",  gp_or,  2);
    gpFunctions[GP_SUB] = new Symbol("sub", gp_sub, 2);
    gpFunctions[GP_XOR] = new Symbol("xor", gp_xor, 2);
    gpFunctions[GP_NOT] = new Symbol("not", gp_not, 1);
    gpFunctions[GP_SR1] = new Symbol("sr1", gp_sr1, 1);
    gpFunctions[GP_SL1] = new Symbol("sl1", gp_sl1, 1);

    ifs.open(argv[1]);
    while(ifs && getline(ifs, buf)){
        if(buf.empty()) continue;
        gpTerminals.push_back(new Symbol(buf));
    }
    ifs.close();
    for(int i=0;i<(int)gpTerminals.size();i++){
        cout << i << " : " << gpTerminals[i]->getSymbol();
        cout << " " << gpTerminals[i]->value.size() << endl;
    }
    
    ifs.open(argv[2]);
    gpAnswer = NULL;
    if(ifs && getline(ifs, buf)){
        gpAnswer = new Symbol(buf);
    }
    ifs.close();

	cout << gpAnswer->getSymbol() << " " << gpAnswer->value.size() << endl;
    
    p.resize(numPopulation);
    for(int i=0;i<(int)p.size();i++){
        p[i] = new Population(numIndividual,
                              numElite,
                              maxDepth,
                              gpFunctions, gpTerminals, gpAnswer);
    }
    // Population(NumIndividual, NumElite, maxDepth);
    /*
	p = new Population(400, 2, 16, gpFunctions, gpTerminals, gpAnswer);
    for(int i=0;i<(int)p->getNumIndividual();i++){
        elite = p->getIndividual(i);
        // cout << elite->getScore() << " : ";
        elite->traverse();
    }
    */
    cout << "Start..." << endl;

    for(int j=0;j<cycleLimit;j++){
        // BinaryTree* getElite(int idx){ return elite[idx]; };
        #ifdef DEBUG
        cout << "==== Cycle: " << j << " ====" << endl;
        cout << " ==== Begin Migration ====" << endl;
        #endif
        
        if(j%migrateInterval == 0){
            random_shuffle(p.begin(), p.end());
            for(int i=1;i<(int)p.size();i++){
                for(int k=0;k<numMigration;k++){
                    p[i]->i_pop_back();
                    p[i]->i_push_back(p[i-1]->getElite(k)->clone());
                }
            }
            for(int k=0;k<numMigration;k++){
                p[0]->i_pop_back();
                p[0]->i_push_back(p[(int)p.size()-1]->getElite(k)->clone());
            }
        }
        #ifdef DEBUG
        cout << " ==== End Migration ====" << endl;
        #endif

        #ifdef DEBUG
        cout << " ==== Begin Presevation ====" << endl;
        #endif
        
        m.clear();
        // #pragma omp parallel for shared(m)
        for(int i=0;i<(int)p.size();i++){
        //p[i]->preserveElite(m, 0);
          p[i]->preserveElite(m, (j/(migrateInterval*10))%2);
          //p[i]->preserveElite(m, 1);
        }
        #ifdef DEBUG
        cout << " ==== End Presevation ====" << endl;
        #endif


        #ifdef DEBUG
        cout << " ==== Begin Selection ====" << endl;
        #endif
        #pragma omp parallel for
        for(int i=0;i<(int)p.size();i++){
            p[i]->selection();
        }
        #ifdef DEBUG
        cout << " ==== End Selection ====" << endl;
        #endif

        #ifdef DEBUG
        cout << " ==== Begin Mutation ====" << endl;
        #endif
        #pragma omp parallel for
        for(int i=0;i<(int)p.size();i++){
            p[i]->mutation();
        }

        #ifdef DEBUG
        cout << " ==== End Mutation ====" << endl;
        #endif

        #ifdef DEBUG
        cout << " ==== Begin CrossOver ====" << endl;
        #endif
        #pragma omp parallel for
        for(int i=0;i<(int)p.size();i++){
            p[i]->crossOver();
        }
        #ifdef DEBUG
        cout << " ==== End CrossOver ====" << endl;
        #endif

        #ifdef DEBUG
        cout << " ==== Begin Elite ====" << endl;
        #endif

        if(j % dispInterval == 0){
            cout << "==== Cycle: " << j << " ====" << endl;
            for(int i=0;i<(int)p.size();i++){
                cout << "P[" << i << "] " << endl;
                for(int k=0;k<numElite;k++){
                    cout << "[" << k << "] S:" << p[i]->getElite(k)->getScore();
                    cout << " N:" << p[i]->getElite(k)->countNode();
                    cout << " ";
                    cout << p[i]->getElite(k)->traverse() << endl;
                }
            }
        }
        #ifdef DEBUG
        cout << " ==== End Elite ====" << endl;
        #endif
    }
    /*
    cout << "score: " << gp_evaluation(gpAnswer->value, tree_output)
         << " node: " << p->getIndividual(0)->getNodeNum()
         << " Depth: " << p->getIndividual(0)->getHeight() << endl;
    */

    return 0;
};

