#include <bits/stdc++.h>

using namespace std;

/**
 * ********************************************************************
 * @brief: 类
 **********************************************************************
 */

class Edge{
public:
    int s;
    int l;
    int id;
    uint64_t channel = 0;

    Edge() {}
    Edge(int s, int l, int id){
        this->id = id;
        this->s = s;
        this->l = l;
    }

    void add(int start, int end){
        uint64_t mask = ((1ULL<<(end - start+1)) -1) << (start-1);
        this->channel |= mask;
    }

    bool ifNULL(int start, int end) const {
        uint64_t mask = ((1ULL << (end - start + 1)) - 1) << (start - 1);
        return (this->channel & mask) == 0;
    }
};

class Node{
public:
    int id;
    unordered_set<int> neighbours;

    Node() {}
    Node(int id){this->id = id;}
};

/**
 * ********************************************************************
 * @brief: 变量
 **********************************************************************
 */

// 超参
const int MAX_EDGE_SIZE = 30;
const string OUTPUT_FILE = R"(..\src\my_case\case1.in)";

int N, M, BusC; // 节点, 边的个数, 业务数量
vector<vector<set<int>>> Graph;
ofstream outFile;
unordered_map<int, Edge> EdgeId_Edge;
unordered_map<int, Node> NodeId_Node;

/**
 * ********************************************************************
 * @brief: 函数
 **********************************************************************
 */

// 写入数据
void writeIntoTxt(const vector<int>& data){
    for(auto i: data) outFile << i << " ";
    outFile << endl;
}

// 生成边的ID, 并随机打乱
void createAndRandEdge(vector<int>& EdgeID){
    for(int i =1; i<=EdgeID.size(); i++) EdgeID[i-1] = i;
    for(int i =0; i<EdgeID.size(); i++){
        int r = rand() % EdgeID.size();
        int t = EdgeID[i];
        EdgeID[i] = EdgeID[r];
        EdgeID[r] = t;
    }
}

void createEdgeInGraph(){
    // 创建Node，并且随机打乱
    vector<int> NodeID(N);
    for(int i =1; i<=NodeID.size(); i++) NodeID[i-1] = i;
    for(int i =0; i<NodeID.size(); i++){
        int r = rand() % NodeID.size();
        int t = NodeID[i];
        NodeID[i] = NodeID[r];
        NodeID[r] = t;
    }

    int p = 0;
    // 保证图联通
    for(int i =0; i<N-1; i++){
        writeIntoTxt(vector<int>{NodeID[i],NodeID[i+1]});
        EdgeId_Edge[i+1] =  Edge(NodeID[i], NodeID[i+1], i+1);
        NodeId_Node[NodeID[i]].neighbours.insert(i+1);
        NodeId_Node[NodeID[i+1]].neighbours.insert(i+1);
        p++;
    }

    while(p<M){
        p++;
        int r1 =rand() % NodeID.size();
        int r2 =rand() % NodeID.size();
        // 无自环
        if(r1==r2) r2 = (r2+1) % NodeID.size();
        EdgeId_Edge[p] =  Edge(NodeID[r1], NodeID[r2], p);
        NodeId_Node[NodeID[r1]].neighbours.insert(p);
        NodeId_Node[NodeID[r2]].neighbours.insert(p);
        writeIntoTxt(vector<int>{NodeID[r1],NodeID[r2]});
    }
}

void printLow40Bits(uint64_t number) {
    const int bitsToPrint = 40;

    uint64_t low40Bits = number & ((1ULL << bitsToPrint) - 1);

    // 将低 40 位转换为二进制字符串
    bitset<bitsToPrint> bitsetLow40(low40Bits);
    string bitString = bitsetLow40.to_string();

    cout << bitString << endl;
}

void createBus(int num){
    // 边的个数
     int edgeSize = rand() % MAX_EDGE_SIZE;

    // 通道起始和结束
    int L = 1 + rand() % (40-1);
    int w = rand()%10 < 6 ? 1 : 1+rand()%(40-L); // 宽度为1的概率较大
    int R = L + w;

    // 输出的边的顺序
    vector<int> outputEdge;

    // 起始节点
    int start;
    int s;

    int count = 0;
    // 随机游走
    while(count < 10000) {
        count++;
        int p = 0;
        // 起始节点
        start = rand() % N + 1;
        s = start;

        set<int> nodes;
        nodes.insert(start);
        while (p < edgeSize) {
            p++;
            bool isOK = false;
            for (int edgeId: NodeId_Node[s].neighbours) {
                Edge &edge = EdgeId_Edge[edgeId];
                if(nodes.count(s == edge.s ? edge.l : edge.s) != 0) continue;
                if (edge.ifNULL(L, R)) {
                    edge.add(L, R);
                    outputEdge.emplace_back(edgeId);
                    s = s == edge.s ? edge.l : edge.s;
                    nodes.insert(s);
                    isOK = true;
                }
                if (isOK) break;
            }
            if (!isOK) break;
        }

        if(outputEdge.empty()) continue;
        else break;
    }

    // 该业务的价值
    int v = rand() % 100000 + 1;
    // int v = R - L + 1;

    writeIntoTxt(vector<int>{start, s, (int)outputEdge.size(), L, R, v});
    writeIntoTxt(outputEdge);
    cout << "finish" << num+1  << endl;
}

int main(){
    outFile.open(OUTPUT_FILE);

    cout << "Input N and M (for example: 10 200)" << endl;
    cin >> N >> M;

    assert(N>=2 && N<=200);
    assert(M>=1 && N<=1000);
    assert(M>=N-1);

    // 写入节点和边
    writeIntoTxt(vector<int>{N, M});

    // 写入变通道次数, 随机生成
    vector<int> P(N);
    for(int i =0; i < N; i++) {
        int r = rand()%10 < 6 ? rand()%2 : rand()%(20+1); // 让0和1的概率大一点
        assert(r>=0 && r<=20);
        P[i] = r;
    }
    writeIntoTxt(P);

    // 生成边的ID
    vector<int> EdgeId(M);
    // 随机打乱
    // CreateAndRandEdge(EdgeId);
    for(int i =1; i<= EdgeId.size(); i++) EdgeId[i-1] = i;

    //创建节点
    for(int i =0 ;i< N;i++) NodeId_Node[i+1] = Node(i+1);

    // 在图上添加边
    createEdgeInGraph();

    // 输入业务数量
    cout << "Input business nums (for example: 20)" << endl;
    cin >> BusC;
    assert(BusC>=0 && BusC<=5000);
    writeIntoTxt(vector<int>{BusC});

    // 生成业务
    for(int i=0;i < BusC; i++) createBus(i);

    // 打印查看稀疏度
    for (auto& p : EdgeId_Edge) {
        Edge& e = p.second;
        printLow40Bits(e.channel);
    }

    // 生成0场景
    writeIntoTxt(vector<int>{0});

    return 0;
}