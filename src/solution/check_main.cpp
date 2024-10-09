#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cstdint>

using namespace std;

using vi = vector<int>;

/**
 * ********************************************************************
 * @brief: 类
 **********************************************************************
 */

class Node {
public:
    int id;
    int p;
    unordered_set<int> neighbors;

    int pCopy;

    Node() {};
    Node(int id, int p): id(id), p(p), pCopy(p) {};

    void back() {
        p = pCopy;
    }
};

class Edge {
public:
    int id;
    int u, v, capacity;
    uint64_t channels0 = 0;  // 记录当前占用，从下标0开始
    uint64_t channels1 = 0;  // 记录重复占用，从下标0开始
    bool isDown = false;

    uint64_t channelsCopy = 0;

    Edge() {};
    Edge( int id, int u, int v, int capacity): id(id), u(u), v(v), capacity(capacity) {};

    void back() {
        isDown = false;
        channels0 = channelsCopy;
        channels1 = 0;
    };

    void costChannels(int start, int width) {
        uint64_t mask = ((1ULL << width) - 1) << (start - 1);
        if ((mask & channels0) != 0)
            channels1 |= (mask & channels0);
        channels0 |= mask;
    }

    void releaseChannels(int start, int width) {
        uint64_t mask = ((1ULL << width) - 1) << (start - 1);
        channels0 &= ~mask;
        channels0 |= channels1;
        channels1 &= (mask ^ channels1);
    }
};

class Business {
public:
    int start, end, width, value, id;
    vector<int> edgePath;
    vector<int> nodePath;
    vector<int> costStartChannels;
    unordered_set<int> usedPNode;

    vector<int> nodePathCopy;
    vector<int> edgePathCopy;
    vector<int> costStartChannelsCopy;

    Business() {};
    Business(int id, int start, int end, int sW, int eW, int value, int lenEdge):
            id(id), start(start), end(end), width(eW-sW+1), value(value) {

        for (int i = 0; i < lenEdge; ++i) {
            costStartChannels.emplace_back(sW);
        }
        costStartChannelsCopy = costStartChannels;
    }

    void back() {
        costStartChannels = costStartChannelsCopy;
        edgePath = edgePathCopy;
        nodePath = nodePathCopy;
        usedPNode.clear();
    }

    void updateP() {
        usedPNode.clear();
        for (int i = 1; i < costStartChannels.size(); ++i) {
            if (costStartChannels[i-1] != costStartChannels[i])
                usedPNode.insert(nodePath[i]);
        }
    }
};

/**
 * ********************************************************************
 * @brief: 变量
 **********************************************************************
 */
// 测试样例个数
int caseNum = 1;
const string CASE_FILE_PREFIX = R"(..\src\semifinal_case\case)";
const string CASE_FILE_SUFFIX = ".in";
const string OUT_FILE_PREFIX = R"(..\src\output\out)";
const string OUT_FILE_SUFFIX = ".txt";
const int capacity = 40;
const float MAX_SIMILARITY = 0.4999;

float score;
int N, M;
int businessNum;
ifstream inFile;
ofstream outFile;
int currentValue = 0;
int ORIGIN_VALUE = 0;
int fixedCount = 0;
int fixedValue = 0;
unordered_map<int, Node> NodeID_Node;
unordered_map<int, Edge> EdgeID_Edge;
unordered_map<int, Business> BusID_Bus;
vector<vector<unordered_set<int>>> Graph;
unordered_set<int> badBusIds;
vector<vector<int>> scenes; // 场景

/**
 * ********************************************************************
 * @brief: 一般函数
 **********************************************************************
 */

int findElementIndex(const vector<int>& vec, const int& element) {
    auto it = find(vec.begin(), vec.end(), element);
    if (it != vec.end()) {
        return distance(vec.begin(), it);
    }
    else {
        return -1;
    }
}

void setNodePath(Business& bus) {
    int nodeId = bus.start;
    assert(findElementIndex(bus.nodePath, nodeId) == -1);
    bus.nodePath.emplace_back(nodeId);
    bus.nodePathCopy.emplace_back(nodeId);

    for (int edgeId : bus.edgePath) {
        Edge& edge = EdgeID_Edge[edgeId];
        int u = edge.u;
        int v = edge.v;
        assert(nodeId == u or nodeId == v);
        if (nodeId == u) {
            bus.nodePath.emplace_back(v);
            bus.nodePathCopy.emplace_back(v);
            nodeId = v;
        } else {
            bus.nodePath.emplace_back(u);
            bus.nodePathCopy.emplace_back(u);
            nodeId = u;
        }
    }
}

void costChannels(const Business& bus) {
    for (int i = 0; i < bus.edgePath.size(); ++i) {
        int edgeId = bus.edgePath[i];
        Edge& edge = EdgeID_Edge[edgeId];

        uint64_t mask = ((1ULL << bus.width) - 1) << (bus.costStartChannels[i] - 1);
        assert((mask & edge.channels0) == 0);

        edge.costChannels(bus.costStartChannels[i], bus.width);
        edge.channelsCopy = edge.channels0;
    }
}

void back() {
    for (auto& entry : BusID_Bus) {
        Business& bus = entry.second;
        bus.back();
    }
    for (auto& entry : NodeID_Node) {
        Node& node = entry.second;
        node.back();
    }
    for (auto& entry : EdgeID_Edge) {
        Edge& edge = entry.second;
        edge.back();
    }

    badBusIds.clear();

    score += float(currentValue) / float(ORIGIN_VALUE) * 10000;
    currentValue = ORIGIN_VALUE;
}

/**
 * Jaccard Similarity: 计算两个序列的相似度（交集除以并集）
 */
float JaccardSimilarity(const vector<int>& seq1, const vector<int>& seq2) {
    if (seq1.empty() || seq2.empty()) return 0;

    unordered_set<int> set1(seq1.begin(), seq1.end());
    unordered_set<int> set2(seq2.begin(), seq2.end());
    // 计算交集
    size_t intersectionSize = 0;
    for (const auto& elem : set1) {
        if (set2.find(elem) != set2.end()) {
            ++intersectionSize;
        }
    }
    // 计算并集
    size_t unionSize = set1.size() + set2.size() - intersectionSize;
    // 返回交集大小除以并集大小（类型转换为浮点数）
    return float(intersectionSize) / float(unionSize);
}

/**
 * ********************************************************************
 * @brief: 数据读写 + 数据检查
 **********************************************************************
 */

//void out(const vector<int>& data, bool lineBreak=true, bool onlyFile=false) {
//    if (onlyFile) {
//        for (const int& n : data) outFile << n << " ";
//        if (lineBreak) outFile << endl;
//    } else {
//        for (const int& n : data) {
//            outFile << n << " ";
//            cout << n << " ";
//        }
//
//        if (lineBreak) {
//            outFile << endl;
//            cout << endl;
//        }
//    }
//}

void out(const vector<int>& data, bool lineBreak=true) {
    for (const int& n : data) {
        outFile << n << " ";
        cout << n << " ";
    }

    if (lineBreak) {
        outFile << endl;
        cout << endl;
    }
}

void loadGraph() {
    inFile >> N >> M;
    assert(N >= 2 and N <= 200);
    assert(M >= 1 and M <= 1000);
    out(vi{N, M});

    // nodes
    vi ps;
    for (int i = 1; i <= N; ++i) {
        int p;

        inFile >> p;
        assert(p >= 0 and p <= 20);
        ps.emplace_back(p);

        Node node(i, p);
        NodeID_Node.emplace(i, node);
    }
    out(ps);

    // edges
    Graph.resize(N + 1, vector<unordered_set<int>>(N + 1));
    for (int i = 1; i <= M; ++i) {
        int u, v;

        inFile >> u >> v;
        assert(u >= 1 and u <= N);
        assert(v >= 1 and v <= N);
        assert(u != v);
        out(vi{u, v});

        Edge edge(i, u, v, capacity);
        EdgeID_Edge.emplace(i, edge);

        Graph[u][v].emplace(i);
        Graph[v][u].emplace(i);

        NodeID_Node[u].neighbors.emplace(v);
        NodeID_Node[v].neighbors.emplace(u);
    }

    // businesses
    inFile >> businessNum;
    assert(businessNum >= 1 and businessNum <= 5000);
    out(vi{businessNum});

    for (int i = 1; i <= businessNum; ++i) {
        int start, end, lenEdge, sW, eW,  v;

        inFile >> start >> end >> lenEdge >> sW >> eW >> v;
        assert(1 <= start and start <= N);
        assert(1 <= end and end <= N);
        assert(1 <= sW and sW <= 40);
        assert(1 <= eW and eW <= 40);
        assert(sW <= eW);
        assert(0 <= v and v <= 100000);
        out(vi{start, end, lenEdge, sW, eW, v});

        Business bus(i, start, end, sW, eW, v, lenEdge);
        vi edgeIds;
        for (int j = 0; j < lenEdge; ++j) {
            int edgeId;

            inFile >> edgeId;
            assert(1 <= edgeId and edgeId <= M);
            edgeIds.emplace_back(edgeId);

            assert(findElementIndex(bus.edgePath, edgeId) == -1);
            bus.edgePath.emplace_back(edgeId);
            bus.edgePathCopy.emplace_back(edgeId);
        }
        out(edgeIds);

        ORIGIN_VALUE += bus.value;
        setNodePath(bus);
        costChannels(bus);
        BusID_Bus.emplace(i, bus);
    }
    currentValue = ORIGIN_VALUE;
}

/**
 * @brief: 加载原有的断边
 */
void loadOriginCase(vector<vector<int>>& scenes) {
    int sceneNum;
    inFile >> sceneNum;

    int n = int(scenes.size());
    for (int i = n; i < sceneNum + n; ++i) {
        int badEdgeNum;
        inFile >> badEdgeNum;

        vector<int> scene;
        for (int _ = 0; _ < badEdgeNum; ++_) {
            int badEdgeId;
            inFile >> badEdgeId;
            assert(findElementIndex(scene, badEdgeId) == -1); // 输出的测试序列中断边重复
            assert(badEdgeId >= 1 and badEdgeId <= M);
            scene.emplace_back(badEdgeId);
        }
        assert(scene.size() == badEdgeNum);
        scene.emplace_back(-1);
        scenes.emplace_back(scene);
    }
}

/**
 * @brief: 加载自定义的断边
 */
void loadCustomCase(vector<vector<int>>& scenes) {
    int sceneNum;
    cin >> sceneNum;

    int n = int(scenes.size());
    for (int i = n; i < sceneNum + n; ++i) {
        int badEdgeNum;
        cin >> badEdgeNum;

        vector<int> scene;
        for (int _ = 0; _ < badEdgeNum; ++_) {
            int badEdgeId;
            cin >> badEdgeId;
            assert(findElementIndex(scene, badEdgeId) == -1); // 输出的测试序列中断边重复
            assert(badEdgeId >= 1 and badEdgeId <= M);
            scene.emplace_back(badEdgeId);
        }

        // 检查相似度
        bool valid = all_of(scenes.begin(), scenes.end(), [&](vector<int>& other) {
            return JaccardSimilarity(scene, other) <= MAX_SIMILARITY;
        });
        assert(valid);

        assert(scene.size() == badEdgeNum);
        scene.emplace_back(-1);
        scenes.emplace_back(scene);
    }
}

void loadCase() {
    loadCustomCase(scenes);
    loadOriginCase(scenes);
}

void loadFixedResult(const vector<int>& scene, vector<vector<vector<vi>>>& allSceneResults) {
    int badEdgeNum = int(scene.size()) - 1;
    vector<vector<vi>> oneSceneResult(badEdgeNum);

    // 循环断边
    for (int en = 0; en < badEdgeNum; ++en) {
        int fixedCount, busId, ePathSize;

        cin >> fixedCount;

        vector<vi> oneEdgeResult(fixedCount); // 一条边坏了的结果
        // 循环一条断边所修复好的所有业务
        for (int i = 0; i < fixedCount; ++i) {
            cin >> busId >> ePathSize;

            vi path = {busId}; // 第一个用来记录busId
            // 循环每个业务的新路径
            for (int j = 0; j < ePathSize; ++j) {
                int edgeId, startChannel, endChannel;

                cin >> edgeId >> startChannel >> endChannel;

                path.emplace_back(edgeId);
                path.emplace_back(startChannel);
                path.emplace_back(endChannel);
            }
            assert(ePathSize == (path.size() - 1) / 3);
            oneEdgeResult[i] = path;
        }
        oneSceneResult[en] = oneEdgeResult;
    }
    allSceneResults.emplace_back(oneSceneResult);
}

void updateBusiness(Business& bus, const vector<int>& newPath) {
    int edgePathSize = int(newPath.size() - 1) / 3;
    int start = bus.start;

    vector<int> newNodePath{ start } ;
    vector<int> newEdgePath;
    vector<int> newCostStartChannels;
    unordered_set<int> visitedNodes{ start };

    for (int i = 0; i < edgePathSize; ++i) {
        int eId = newPath[i * 3 + 1];
        Edge& e = EdgeID_Edge[eId];
        assert(eId >= 1 and eId <= M); // 路径的边ID不对
        assert(!e.isDown); // 路径经过了已经中断的光纤
        assert(e.u == start or e.v == start); // 路径不通

        start = start == e.u ? e.v : e.u;
        assert(visitedNodes.find(start) == visitedNodes.end()); // 路径成环
        visitedNodes.emplace(start);
        if (i == edgePathSize - 1) assert(start == bus.end);  // 路径不通

        int startChannel = newPath[i * 3 + 2];
        int endChannel = newPath[i * 3 + 3];
        assert(bus.width == endChannel - startChannel + 1); // 路径上的业务宽度不一致
        assert(startChannel <= endChannel); // 业务所使用的通道ID不对
        assert(startChannel >= 1 and startChannel <= capacity); // 业务所使用的通道ID不对
        assert(endChannel >= 1 and endChannel <= capacity); // 业务所使用的通道ID不对

        newCostStartChannels.emplace_back(startChannel);
        newEdgePath.emplace_back(eId);
        newNodePath.emplace_back(start);
    }

    bus.edgePath = newEdgePath;
    bus.nodePath = newNodePath;
    bus.costStartChannels = newCostStartChannels;
    bus.updateP();
}

vector<int> findFailedBusIds(const int& failedEdgeId, const unordered_set<int>& abandonedBusIds) {
    vector<int> failedBusIds;
    for (auto& entry : BusID_Bus) {
        Business& bus = entry.second;
        if (abandonedBusIds.find(bus.id) == abandonedBusIds.end() && find(bus.edgePath.begin(), bus.edgePath.end(), failedEdgeId) != bus.edgePath.end()) {
            failedBusIds.emplace_back(bus.id);
            currentValue -= bus.value;
        }
    }
    return failedBusIds;
}

void releaseResource(const Business& newBus, const Business& oldBus) {
    for (int i = 0; i < oldBus.edgePath.size(); ++i) {
        int edgeId = oldBus.edgePath[i];
        Edge& edge = EdgeID_Edge[edgeId];
        edge.releaseChannels(oldBus.costStartChannels[i], oldBus.width);
        if (i > 0 && oldBus.costStartChannels[i-1] != oldBus.costStartChannels[i]) {
            int nodeId = oldBus.nodePath[i];
            if (newBus.usedPNode.find(nodeId) == newBus.usedPNode.end()) {
                Node& node = NodeID_Node[nodeId];
                ++node.p;
            }
        }
    }
}

void updateGraph(const Business& newBus, const Business& oldBus) {
    for (int i = 0; i < newBus.edgePath.size(); ++i) {
        int edgeId = newBus.edgePath[i];
        Edge& edge = EdgeID_Edge[edgeId];

        // check
        uint64_t mask = ((1Ull << newBus.width) - 1) << (newBus.costStartChannels[i] - 1);
        uint64_t channelsCopy = edge.channels0;
        int idx = findElementIndex(oldBus.edgePath, edgeId);
        if (idx != -1) {
            uint64_t _mask =  ((1ULL << oldBus.width) - 1) << (oldBus.costStartChannels[idx] - 1);
            channelsCopy &= ~_mask;
        }
        assert((channelsCopy & mask) == 0); // 业务所使用的通道被占用

        edge.costChannels(newBus.costStartChannels[i], newBus.width);
        if (i > 0 && newBus.costStartChannels[i-1] != newBus.costStartChannels[i]) {
            int nodeId = newBus.nodePath[i];
            if (oldBus.usedPNode.find(nodeId) == oldBus.usedPNode.end()) {
                Node& node = NodeID_Node[nodeId];
                --node.p;
                assert(node.p >= 0);  // 节点变通道次数不够
            }
        }
    }
}

void handleBadEdge(const int& edgeId, const vector<vi>& oneEdgeResult) {
    assert(edgeId >= 1 and edgeId <= M); // 路径的边ID不对
    assert(EdgeID_Edge.find(edgeId) != EdgeID_Edge.end()); // 路径的边ID不对

    unordered_set<int> visitedBusIds;
    vector<int> failedBusIds = findFailedBusIds(edgeId, badBusIds);

    vector<Business> releasedBusinesses;
    vector<Business> fixedBusinesses;

    // set isDown
    Edge& edge = EdgeID_Edge[edgeId];
    edge.isDown = true;

    for (const vector<int>& newPath : oneEdgeResult) {
        assert((newPath.size() - 1) % 3 == 0);

        int busId = newPath[0];
        assert(busId >= 1 && busId <= businessNum); // 输出的业务ID不对
        assert(visitedBusIds.find(busId) == visitedBusIds.end());  // 单次输出的业务ID重复
        assert(badBusIds.find(busId) == badBusIds.end()); // 输出已经废弃的业务
        assert(findElementIndex(failedBusIds, busId) != -1); // 输出未受到此次光纤中断影响的业务ID
        visitedBusIds.emplace(busId);

        Business& bus = BusID_Bus[busId];
        Business oldBusCopy = bus;  // 复制

        updateBusiness(bus, newPath);
        updateGraph(bus, oldBusCopy);

        fixedBusinesses.emplace_back(bus);
        releasedBusinesses.emplace_back(oldBusCopy);
    }

    for (int busId : failedBusIds) {
        if (visitedBusIds.find(busId) == visitedBusIds.end()) badBusIds.emplace(busId);
    }

    for (int i = 0; i < releasedBusinesses.size(); ++i) {
        releaseResource(fixedBusinesses[i], releasedBusinesses[i]);
        currentValue += fixedBusinesses[i].value;
        fixedValue += fixedBusinesses[i].value;
        fixedCount += 1;
    }
}

/**
 * 将算发输出写入文件
 */
void writeAllSceneResults(const vector<vector<vector<vi>>>& allSceneResults) {
    for (auto& sceneResult : allSceneResults) {
        for (auto& oneEdgeResult : sceneResult) {
            outFile << oneEdgeResult.size() << endl;
            for (auto& path : oneEdgeResult) {
                assert((path.size() - 1) % 3 == 0);
                int edgeNum = int((path.size() - 1) / 3);
                outFile << path[0] << " " << edgeNum << endl;
                for (int i = 1; i < path.size(); ++i) {
                    outFile << path[i] << " ";
                }
                outFile << endl;
            }
        }
    }
}

/**
 * 验证算法输出的合理性
 */
void validResults(const vector<vector<int>>& scenes, const vector<vector<vector<vi>>>& allSceneResults) {
    assert(scenes.size() == allSceneResults.size());
    for (int sn = 0; sn < scenes.size(); ++sn) {

        const vector<int>& scene = scenes[sn];
        const vector<vector<vi>>& sceneResult = allSceneResults[sn];
        assert(scene.size() - 1 == sceneResult.size());

        for (int en = 0; en < scene.size(); ++en) {
            if (en == scene.size() - 1) {
                assert(scene[en] == -1);
                back(); // 重置
            } else {
                handleBadEdge(scene[en], sceneResult[en]);
            }
        }
    }
}

/**
 * ********************************************************************
 * @brief: 主函数
 **********************************************************************
 */

void runSolution() {
    loadGraph();

    // 加载自定义的场景
    double startTimeG = clock();
    loadCustomCase(scenes);
    double endTimeG = clock();
    // 加载官方场景
    loadOriginCase(scenes);

    // run solution
    double startTimeH = clock();
    vector<vector<vector<vi>>> results;  // {每一场景{每一边{每一业务{每一修复路径}}}}
    out(vi{int(scenes.size())});
    for (const vector<int>& scene : scenes) {
        out(scene);
        loadFixedResult(scene, results);
    }
    double endTimeH = clock();

    cout << "running time of generation: " << (endTimeG - startTimeG) / 1000 << "s" << endl;
    cout << "running time of handle: " << (endTimeH - startTimeH) / 1000 << "s" << endl;


    // write to file
    writeAllSceneResults(results);

    // 验证输出的有效性
    validResults(scenes, results);

    cout << "score: " << score << endl;
    cout << "fixedCount: " << fixedCount << endl;
    cout << "fixedValue: " << fixedValue << endl;
}

int main(int argc, char **argv) {
    // 读取参数
    if (argc > 1) {
        caseNum = stoi(argv[1]);  // 将字符串转换为整数
    }

    // 初始化文件流
    inFile.open(CASE_FILE_PREFIX + to_string(caseNum) + CASE_FILE_SUFFIX);
    outFile.open(OUT_FILE_PREFIX + to_string(caseNum) + OUT_FILE_SUFFIX);

    runSolution();

    inFile.close();
}
