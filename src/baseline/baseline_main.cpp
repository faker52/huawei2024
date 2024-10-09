//#define NDEBUG // 发布模式，禁用assert

#include <iostream>
#include <set>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <deque>
#include <memory>
#include <limits>
#include <bitset>
#include <queue>

using namespace std;

/**
* *******************************************************************
* @brief: 类
* *******************************************************************
*/

class Node {
public:
    int ID;
    int leftChangeTime;
    vector<int> neighborEIds;
};

class Edge {
public:
    int ID;
    int uNId, vNId;
    bool isDown = false;
    uint64_t channels = 0;

    void costChannels(int start, int width) {
        uint64_t mask = ((1ULL << width) - 1) << (start - 1);
        assert((channels & mask) == 0);
        channels |= mask;
    }

    int otherNId(int nId) const {
        return nId == uNId ? vNId : uNId;
    }
};

class Path {
public:
    int startNId;
    int width;

    vector<int> nPath;
    vector<int> ePath;
    vector<int> startC;

    unordered_set<int> usePNodes; // 使用变节点能力的节点

    bool hasEdge(int eId) const {
        return (find(ePath.begin(), ePath.end(), eId) != ePath.end());
    }

    uint64_t getMask(int idx) const {
        return ((1ULL << width) - 1) << (startC[idx] - 1);
    }

    bool empty() const {
        return ePath.empty() || startC.empty();
    }

    void calUsePNodes() {
        for (int i = 0; i < ePath.size() - 1; ++i) {
            if (startC[i] != startC[i+1])
                usePNodes.emplace(nPath[i+1]);
        }
    }
};

class Service {
public:
    int ID;
    int startNId, endNId;
    int width;
    int value;
    bool isAlive = true;
    Path path;

    void set(int ID, int sNId, int eNId, int leftC, int rightC, int value) {
        this->ID = ID;
        this->startNId = sNId;
        this->endNId = eNId;
        this->width = rightC - leftC + 1;
        this->value = value;

        path.width = this->width;
        path.startNId = startNId;
    }
};

class ANode {
public:
    int idx = 1; // 位于路径中的位置，从1开始
    int nId;
    int fromEId;
    uint64_t fromC; // 边通道
    shared_ptr<ANode> lastANode;

    ANode(int nId, int fromEId = -1, uint64_t fromC = 0, const shared_ptr<ANode>& lastANode = nullptr)
            : nId(nId), fromEId(fromEId), fromC(fromC), lastANode(lastANode) {
        if (lastANode != nullptr) idx = lastANode->idx + 1;
    };

    bool isNodeInPath(int nodeId) const {
        const ANode* ptr = lastANode.get();
        while (ptr != nullptr) {
            if (ptr->nId == nodeId) return true;
            ptr = ptr->lastANode.get();
        }
        return false;
    }
};

/**
* *******************************************************************
* @brief: 变量
* *******************************************************************
*/
// 变量
int nodeNum;
int edgeNum;
int serviceNum;

int totalValue = 0; // 总价值
int currentValue = 0; // 当前价值
vector<Node> originNodes;
vector<Edge> originEdges;
vector<Service> originServices;
unordered_map<int, vector<shared_ptr<ANode>>> shortestPaths; // 记录全局最短路径

// 常量
const unsigned seed = 0;
const int INF = 1e9;
const int capacity = 40;
const float PATH_LEN_EXPANSION = 1.2;

/**
* *******************************************************************
* @brief: 函数
* *******************************************************************
*/

void printLow40Bits(uint64_t number) {
    const int bitsToPrint = 40;

    uint64_t low40Bits = number & ((1ULL << bitsToPrint) - 1);

    // 将低 40 位转换为二进制字符串
    bitset<bitsToPrint> bitsetLow40(low40Bits);
    string bitString = bitsetLow40.to_string();

    cout << bitString << endl;
}

void loadEnv() {
    cin >> nodeNum >> edgeNum;

    // 获取每个节点变通道数
    originNodes.resize(nodeNum + 1);
    for (int i = 1; i <= nodeNum; ++i) {
        cin >> originNodes[i].leftChangeTime;
        originNodes[i].leftChangeTime = 0; // baselines
        originNodes[i].ID = i;
    }

    // 获取边
    originEdges.resize(edgeNum + 1);
    for (int i = 1; i <= edgeNum; ++i) {
        int uNId, vNId;
        cin >> uNId >> vNId;

        originNodes[uNId].neighborEIds.emplace_back(i);
        originNodes[vNId].neighborEIds.emplace_back(i);

        originEdges[i].ID = i;
        originEdges[i].uNId = uNId;
        originEdges[i].vNId = vNId;
    }

    // 获取业务
    cin >> serviceNum;
    originServices.resize(serviceNum + 1);
    for (int i = 1; i <= serviceNum; ++i) {
        int sNId, eNId, len, leftC, rightC, value;
        cin >> sNId >> eNId >> len >> leftC >> rightC >> value;

        Service& service = originServices[i];
        service.set(i, sNId, eNId, leftC, rightC, value);

        int srcNId = sNId;
        service.path.nPath.emplace_back(srcNId);

        // 初始化路径
        for (int _ = 0; _ < len; ++_) {
            int edgeId;
            cin >> edgeId;

            service.path.ePath.emplace_back(edgeId);
            service.path.startC.emplace_back(leftC);

            srcNId = originEdges[edgeId].otherNId(srcNId);
            service.path.nPath.emplace_back(srcNId);

            originEdges[edgeId].costChannels(leftC, service.width);
        }

        // 计算价值
        totalValue += value;
    }
}

int findElementIndex(const vector<int>& vec, int element) {
    auto it = find(vec.begin(), vec.end(), element);
    if (it != vec.end()) {
        return distance(vec.begin(), it);
    }
    else {
        return -1;
    }
}

int findContinueIndex(uint64_t pre, int w) {
    uint64_t mask = (1ULL << w) - 1;
    for (int i = 0; i <= 40 - w; ++i) {
        if ((pre & mask) == 0) {
            return i + 1;
        }
        pre >>= 1;
    }
    return -1;
}

uint64_t canConnect(uint64_t lastChannels, uint64_t nextChannels, int width, bool isUsePNode, int leftChangeTime) {
    assert(leftChangeTime >= 0);
    // 不可变通道
    if (leftChangeTime == 0 && !isUsePNode)
        nextChannels |= lastChannels;
    if (findContinueIndex(nextChannels, width) == -1) return numeric_limits<uint64_t>::max(); // 前后通道不通
    else return nextChannels;
}

vector<int> findAllContinueIndex(uint64_t channels, int width)
{
    const uint64_t mask = (1ULL << width) - 1;
    vector<int> startChannels;

    for (int i = 1; i <= 40 - width + 1; ++i)
    {
        if ((channels & mask) == 0)
            startChannels.emplace_back(i);
        channels >>= 1;
    }
    return startChannels;
}

/**
* *******************************************************************
* @brief: tricks
* *******************************************************************
*/

// 根据占用大小排序，越大越前
void orderEdges1(vector<int>& edgeIds, const vector<Edge>& edges, const vector<Node>& nodes) {
    sort(edgeIds.begin(), edgeIds.end(), [&](int a, int b) {
        int aNumOnes = __builtin_popcountll(edges[a].channels);
        int bNumOnes = __builtin_popcountll(edges[b].channels);
        return aNumOnes < bNumOnes;
    });
}

// 根据两端节点的变通道能力总数排序，越小越前
void orderEdges2(vector<int>& edgeIds, const vector<Edge>& edges, const vector<Node>& nodes) {
    sort(edgeIds.begin(), edgeIds.end(), [&](int a, int b) {
        int aChangeTimes = nodes[edges[a].uNId].leftChangeTime + nodes[edges[a].vNId].leftChangeTime;
        int bChangeTimes = nodes[edges[b].uNId].leftChangeTime + nodes[edges[b].vNId].leftChangeTime;
        return aChangeTimes > bChangeTimes;
    });
}

void tricks(vector<Edge>& edges, vector<Node>& nodes) {
    for (Node& node : nodes) {
        sort(node.neighborEIds.begin(), node.neighborEIds.end(), [&](int a, int b) {
            int aNumOnes = __builtin_popcountll(edges[a].channels);
            int bNumOnes = __builtin_popcountll(edges[b].channels);

            int aChangeTimes = nodes[edges[a].uNId].leftChangeTime + nodes[edges[a].vNId].leftChangeTime;
            int bChangeTimes = nodes[edges[b].uNId].leftChangeTime + nodes[edges[b].vNId].leftChangeTime;

            return aNumOnes < bNumOnes;
            if (aChangeTimes == bChangeTimes) return aNumOnes > bNumOnes;
            return aChangeTimes < bChangeTimes;
        });
    }
}

/**
* *******************************************************************
* @brief: 寻路通用函数
* *******************************************************************
*/

void buildPath(const shared_ptr<ANode>& endANode, Path& path) {
    const ANode* ptr = endANode.get();
    while (ptr->fromEId != -1) {
        path.nPath.emplace_back(ptr->nId);
        path.ePath.emplace_back(ptr->fromEId);
        ptr = ptr->lastANode.get();
    }

    // 反转边路径
    reverse(path.ePath.begin(), path.ePath.end());

    // 反转节点路径
    path.nPath.emplace_back(path.startNId);
    reverse(path.nPath.begin(), path.nPath.end());
}

vector<int> getDistance(const Service& serv, const vector<Node>& nodes, const vector<Edge>& edges) {
    unordered_set<int> visited = { serv.endNId };
    vector<int> dis(nodes.size() + 1, INF); // 存储各节点距离终点的距离
    dis[serv.endNId] = 0;

    queue<int> q;
    q.emplace(serv.endNId);

    while (!q.empty()) {
        int src = q.front();
        q.pop();

        for (int eId : nodes[src].neighborEIds) {
            // 跳过坏边
            const Edge& edge = edges[eId];
            if (edge.isDown) continue;

            // 跳过已经访问的节点
            int neighbor = edges[eId].otherNId(src);
            if (visited.find(neighbor) != visited.end()) continue;

            // 跳过部分不通的边
            if (findContinueIndex(edge.channels, serv.width) == -1) continue;

            visited.insert(neighbor);
            dis[neighbor] = dis[src] + 1;
            q.emplace(neighbor);
        }
    }
    return dis;
}

/**
* *******************************************************************
* @brief: 寻路baseline
* 先寻找全局最短路径（考虑通道连接性）
* *******************************************************************
*/

struct Cash {
    int count = 800;
    unordered_map<int, vector<shared_ptr<ANode>>> visitedForward;
    unordered_map<int, vector<shared_ptr<ANode>>> visitedBackward;
};

// 找变通道最少的通道
void buildChannelsBase(uint64_t channels, Path& path) {
    // 找到起始通道编号
    int startChannel = findContinueIndex(channels, path.width);
    assert(startChannel != -1);

    // 重构完整通道路径
    for (int i = 0; i < path.ePath.size(); ++i) {
        path.startC.emplace_back(startChannel);
    }
}

void buildPathBase(const shared_ptr<ANode>& endANodeF, const shared_ptr<ANode>& endANodeB, Path& path) {
    // 添加前向
    const ANode* ptrF = endANodeF.get();
    while (ptrF->fromEId != -1) {
        path.nPath.emplace_back(ptrF->nId);
        path.ePath.emplace_back(ptrF->fromEId);
        ptrF = ptrF->lastANode.get();
    }

    // 反转边路径
    reverse(path.ePath.begin(), path.ePath.end());

    // 反转节点路径
    path.nPath.emplace_back(path.startNId);
    reverse(path.nPath.begin(), path.nPath.end());

    // 添加后向节点
    const ANode* ptrB = endANodeB.get();
    while (ptrB->fromEId != -1) {
        path.ePath.emplace_back(ptrB->fromEId);
        ptrB = ptrB->lastANode.get();
        path.nPath.emplace_back(ptrB->nId);
    }
}

void seekPathBase(Service& serv, vector<Edge>& edges, vector<Node>& nodes, Path& path,
                  deque<shared_ptr<ANode>>& q, Cash& cash, bool isForward) {

    auto& v1 = isForward ? cash.visitedForward : cash.visitedBackward;
    auto& v2 = isForward ? cash.visitedBackward : cash.visitedForward;
    v1.clear();

    const size_t s = q.size();
    for (int i = 0; i < s; ++i) {
//        // 超出循环次数限制，则退出
//        if (cash.count <= 0) return;
//        else cash.count = cash.count - 1;

        shared_ptr<ANode> aNode = q.front();
        q.pop_front();

        const Node& node = nodes[aNode->nId];

        // 检查前后向路径是否可以拼接
        if (v2.find(aNode->nId) != v2.end()) {
            for (const shared_ptr<ANode>& other : v2[aNode->nId]) {
                uint64_t channels = canConnect(aNode->fromC, other->fromC, serv.width, false, 0); // 只允许不变通道
                if (channels != numeric_limits<uint64_t>::max()) { // 路径可通
                    if (isForward) buildPathBase(aNode, other, path);
                    else buildPathBase(other, aNode, path);
                    buildChannelsBase(channels, path);
                    return;
                }
            }
        }

        for (int eId : node.neighborEIds) {
            const Edge& edge = edges[eId];

            // 跳过坏边
            if (edge.isDown) continue;

            // 跳过在路径中节点
            int nextNId = edge.otherNId(aNode->nId);
            if (aNode->isNodeInPath(nextNId)) continue;

            // 检查前后边的通道是否通
            uint64_t channels = canConnect(aNode->fromC, edge.channels, serv.width, false, 0); // 只允许不变通道
            if (channels == numeric_limits<uint64_t>::max())  // 不通
                continue;

            shared_ptr<ANode> nextANode = make_shared<ANode>(nextNId, eId, channels, aNode);
            q.emplace_back(nextANode);
            v1[nextNId].emplace_back(nextANode);
        }
    }
}

void rebuildPathBase(Service& serv, vector<Edge>& edges, vector<Node>& nodes, Path& path, int maxIter=10) {
    // 初始化前向队列
    deque<shared_ptr<ANode>> qf;
    qf.emplace_back(make_shared<ANode>(serv.startNId, -1, serv.path.getMask(0), nullptr));

    // 初始化后向队列
    deque<shared_ptr<ANode>> qb;
    qb.emplace_back(make_shared<ANode>(serv.endNId, -1, serv.path.getMask(0), nullptr));

    // 共享数据
    Cash cash;

//    while (true) {
    for (int i = 0; i < maxIter; ++i) {
        seekPathBase(serv, edges, nodes, path, qf, cash, true);
        if (cash.visitedForward.empty() || cash.count <= 0 || !path.empty()) return;
        seekPathBase(serv, edges, nodes, path, qb, cash, false);
        if (cash.visitedBackward.empty() || cash.count <= 0 || !path.empty()) return;
    }
}


/**
* *******************************************************************
* @brief: 输出
* *******************************************************************
*/

void outputResult(const Service& serv) {
    cout << serv.ID << " " << serv.path.ePath.size() << endl;

    for (int i = 0; i < serv.path.ePath.size(); ++i) {
        int startChannel = serv.path.startC[i];
        cout << serv.path.ePath[i] << " " << startChannel << " " << startChannel + serv.width - 1 << " ";
    }

    cout << endl;
    fflush(stdout);
}

/**
* *******************************************************************
* @brief: 修复
* *******************************************************************
*/

void releaseRes(const Service& serv, vector<Edge>& edges, vector<Node>& nodes) {
    // 释放通道
    for (int i = 0; i < serv.path.ePath.size(); ++i) {
        Edge& edge = edges[serv.path.ePath[i]];
        edge.channels &= ~serv.path.getMask(i);
    }

    // 释放变节点能力
    for (int nId : serv.path.usePNodes) {
        ++nodes[nId].leftChangeTime;
    }
}

void releaseRes(const Service& oldServ, const Service& newServ, vector<Edge>& edges, vector<Node>& nodes) {
    // 释放旧业务占用，而新业务未占用的通道
    const Path& oldPath = oldServ.path;
    const Path& newPath = newServ.path;
    for (int i = 0; i < oldPath.ePath.size(); ++i) {
        Edge& edge = edges[oldServ.path.ePath[i]];
        int idx = findElementIndex(newPath.ePath, oldPath.ePath[i]);

        if (idx == -1) edge.channels &= ~oldServ.path.getMask(i);
        else edge.channels &= ~(oldPath.getMask(i) & newPath.getMask(idx) ^ oldPath.getMask(i));
    }

    // 释放旧业务占用，而新业务未占用的变通道能力
    for (int nId : oldPath.usePNodes) {
        if (newPath.usePNodes.find(nId) == newPath.usePNodes.end()) ++nodes[nId].leftChangeTime;
    }
}

void occupyRes(const Service& serv, vector<Edge>& edges, vector<Node>& nodes, const Path& path) {
    // 占据旧业务通道
    for (int i = 0; i < serv.path.ePath.size(); ++i) {
        Edge& edge = edges[serv.path.ePath[i]];
        edge.channels |= serv.path.getMask(i);
    }

    // 占据新业务通道
    for (int i = 0; i < path.ePath.size(); ++i) {
        Edge& edge = edges[path.ePath[i]];
        edge.channels |= path.getMask(i);
    }

    // 使用变通道能力
    unordered_set<int> unionSet;
    for (int i : path.usePNodes) unionSet.emplace(i);
    for (int i : serv.path.usePNodes) unionSet.emplace(i);
    for (int i : unionSet) --nodes[i].leftChangeTime;
}

vector<int> handleFailedEdge(int failedEId, vector<Service>& services, vector<Edge>& edges, vector<Node>& nodes) {
    // 设置坏边
    edges[failedEId].isDown = true;

    // 找到受影响的业务
    vector<int> affectedServIds;
    for (Service& serv : services) {
        if (serv.isAlive && serv.path.hasEdge(failedEId)) {
            affectedServIds.emplace_back(serv.ID);
            currentValue -= serv.value;
        }
    }

    // 对业务按照价值排序，越高越前
    sort(affectedServIds.begin(), affectedServIds.end(), [&services](int aId, int bId) {
        return services[aId].value > services[bId].value;
    });

    // 存储修复业务的旧版本，以便后续释放旧资源
    vector<Service> fixedOldServices;
    // 成功修复的业务
    vector<int> fixedServIds;

    // 修复业务
    for (int sId : affectedServIds) {
        Service& serv = services[sId];

        // 业务不存活，则跳过
        if (!serv.isAlive) continue;

        // 释放原有业务的资源
        releaseRes(serv, edges, nodes);

        // tricks
//        tricks(edges, nodes);

        // 重新路径规划路径
        Path newPath;
        newPath.startNId = serv.startNId;
        newPath.width = serv.width;
        rebuildPathBase(serv, edges, nodes, newPath);

        // 重新占用资源
        occupyRes(serv, edges, nodes, newPath);

        if (!newPath.empty()) { // 规划成功
            fixedServIds.emplace_back(serv.ID);
            fixedOldServices.emplace_back(serv); // 复制，以便后续释放旧资源

            // 更新业务
            serv.isAlive = true;
            serv.path = newPath;
        } else { // 规划失败
            // 更新业务
            serv.isAlive = false;
        }
    }

    // 释放旧业务资源
    for (Service& oldServ : fixedOldServices) {
        releaseRes(oldServ, services[oldServ.ID], edges, nodes);
    }

    return fixedServIds;
}

void handleCases() {
    // 获取测试场景
    int caseNum;
    cin >> caseNum;

    for (int _ = 0; _ < caseNum; ++_) {

        vector<Service> services = originServices;
        vector<Edge> edges = originEdges;
        vector<Node> nodes = originNodes;

        currentValue = totalValue; // 重置当前价值

        while (true) {
            int failedEId;
            cin >> failedEId;

            if (failedEId == -1) break;

            // 处理坏边
            const vector<int>& fixedServIds = handleFailedEdge(failedEId, services, edges, nodes);

            // 输出结果
            cout << fixedServIds.size() << endl;
            for (int sId : fixedServIds) {
                assert(services[sId].isAlive);
                outputResult(services[sId]);

                currentValue += services[sId].value;
            }
        }
    }
}
/**
* *******************************************************************
* @brief: 生成数据
* *******************************************************************
*/

void generateCases() {
    cout << 0 << endl;
}

/**
* *******************************************************************
* @brief: 主函数
* *******************************************************************
*/

int main() {
    // 设置随机种子
    srand(seed);

    // 加载环境
    loadEnv();

    // 生成数据
    double startTimeG = clock();
//    generateCases();
    double endTimeG = clock();

    // 记录时间
    double startTimeH = clock();
    handleCases();
    double endTimeH = clock();

    cout << "running time of generation: " << (endTimeG - startTimeG) / 1000 << "s" << endl;
    cout << "running time of handle: " << (endTimeH - startTimeH) / 1000 << "s" << endl;

    return 0;
}