#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <iterator>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<map>
#include<algorithm>
#include<queue>
#include <chrono>

using namespace std;
using namespace std::chrono;

struct Edge
{
    int U;
    int V;
    bool operator==(const Edge& e)const
    {
        return U == e.U && V == e.V;
    }

};

bool edgeOp(Edge e1, Edge e2)
{
    if(e1.U < e2.U)
        return true;
    else if(e1.U == e2.U)
        return e1.V < e2.V;
    else
        return false;

}

class EdgeHashFunction
{
public:
    size_t operator()(const Edge edgeToHash) const noexcept
    {
            size_t hash = edgeToHash.U + 10 * edgeToHash.V;
            return hash;
    }
};

struct structureDegree
{
    int Degree;
    int Vertex;
};

struct structureVertex
{
    int Vertex;
    list<structureVertex>::iterator Itr;
};

bool opDeg(structureDegree s1, structureDegree s2)
{
    return s1.Degree > s2.Degree;
}

enum reductionType {bussReductionMethod, crownReductionMethod};



void delta(vector<int> augPath, vector<Edge> maxMatching, vector<Edge>& newMatching)
{
   vector<Edge> augEdgeSet;
   Edge edge;

   for(int i=0; i<augPath.size()-1; i++){
       if(augPath[i] < augPath[i+1]){
           edge.U = augPath[i];
           edge.V = augPath[i+1];
       }
       else {
           edge.U = augPath[i+1];
           edge.V = augPath[i];
       }
       augEdgeSet.push_back(edge);
   }

   sort(maxMatching.begin(), maxMatching.end(), edgeOp);
   sort(augEdgeSet.begin(), augEdgeSet.end(), edgeOp);

//   cout << "-----------------previous maximum matching--------------"<< endl;
//   for(int i=0; i < maxMatching.size(); i++)
//   {
//        cout << maxMatching[i].U <<"-"<< maxMatching[i].V<< endl;
//   }
//   cout << "----------------augmenting path edges -------------------"<< endl;
//   for(int i=0; i<augEdgeSet.size(); i++)
//   {
//       cout << augEdgeSet[i].U << "-" << augEdgeSet[i].V << endl;
//   }


   int i=0, j=0;
    while(i != maxMatching.size() && j != augEdgeSet.size())
    {
//        cout << "i -> " << i << " j -> " << j << endl;
        if(maxMatching[i] == augEdgeSet[j]){
            i++;
            j++;
//            cout << "if part" << endl;
        }
        else if(edgeOp(maxMatching[i], augEdgeSet[j])){
            newMatching.push_back(maxMatching[i]);
            i++;
//            cout << "else if part" << endl;
        }
        else {
            newMatching.push_back(augEdgeSet[j]);
            j++;
//            cout << "else part" << endl;
        }
    }

    if(i == maxMatching.size()){
        while(j != augEdgeSet.size()){
            newMatching.push_back(augEdgeSet[j]);
            j++;
        }
    }
    else{
        while(i != maxMatching.size()){
            newMatching.push_back(maxMatching[i]);
            i++;
        }
    }


    return;
}



void setMinus(vector<int> setA, vector<int> setB, vector<int>& result)
{
    sort(setA.begin(), setA.end());
    sort(setB.begin(), setB.end());

    int i = 0, j = 0;

    while(i != setA.size() && j != setB.size())
    {
        if(setA[i] == setB[j])
        {
            i++;
            j++;
        }
        else if(setA[i]  < setB[j])
        {
            result.push_back(setA[i]);
            i++;
        }
        else {
            j++;
        }
    }

    if(j == setB.size()){
        while(i != setA.size())
        {
            result.push_back(setA[i]);
            i++;
        }
        return;
    }


}


void setIntersection(vector<int> setA, vector<int> setB, vector<int>& result)
{
    int i=0, j=0;
    sort(setA.begin(), setA.end());
    sort(setB.begin(), setB.end());
    while(i != setA.size() && j != setB.size())
    {
        if(setA[i] == setB[j])
        {
            result.push_back(setA[i]);
            i++;
            j++;
        }
        else if(setA[i]  < setB[j])
            i++;
        else {
            j++;
        }
    }
    return;
}




void makeCombiUtil(vector<vector<int> >& ans,
    vector<int>& tmp, vector<int>& data, int n, int left, int k)
{
    // Pushing this vector to a vector of vector
    if (k == 0) {
        ans.push_back(tmp);
        return;
    }

    // i iterates from left to n. First time
    // left will be 1
    for (int i = left; i < n; ++i)
    {
        tmp.push_back(data[i]);
        makeCombiUtil(ans, tmp, data, n, i + 1, k - 1);

        // Popping out last inserted element
        // from the vector
        tmp.pop_back();
    }
}

vector<vector<int> > makeCombi(int k, vector<int> data)
{
    int n = data.size();
    vector<vector<int> > ans;
    vector<int> tmp;
    makeCombiUtil(ans, tmp, data, n, 0, k);
    return ans;
}

void print(vector<int> data)
{
    for(int i=0; i<data.size(); i++)
        cout << data[i] << " ";
    cout << endl;
    return;
}

void print(vector<structureDegree> data)
{
    for(int i=0; i<data.size(); i++)
    {
        cout << data[i].Vertex << "->" << data[i].Degree << " |";
    }
    cout << endl;
    return;
}

void print(vector<Edge> data)
{
    for(int i=0; i<data.size(); i++)
    {
        cout << data[i].U << "-" << data[i].V <<" | ";
    }

    cout <<"****"<< endl;
}


class Graph
{
public:

    int numberOfVertices;
    int numberOfEdges;
    unordered_set<Edge, EdgeHashFunction> edgeSet;
    vector<list<structureVertex>> adjacencyList;
    unordered_set<int> vertexSet;
    vector<int> vertexCover;
    vector<Edge> maximalMatching;
    vector<Edge> independentEdgeSet;
    vector<structureDegree> degreeVector;
    vector<int> degreeIdxVector;

    Graph(int nv, int ne)
    {
        vector<list<structureVertex>> vec1(nv+1);
        vector<int> vec2(nv+1);
        numberOfVertices = nv;
        numberOfEdges = ne;
        adjacencyList = vec1;
        degreeIdxVector = vec2;
        setVertices();
//        cout << nv << " vertex added" << endl;
    }

    void setVertices()
    {
        for(int i=1; i<=numberOfVertices; i++)
        {
            vertexSet.insert(i);
        }
        // cout << "set vertex" << endl;
        return;
    }


    void addEdge(unsigned u, unsigned v)
    {
        structureVertex sVertex1, sVertex2;

        sVertex1.Vertex = v;
        sVertex2.Vertex = u;

        adjacencyList[u].push_back(sVertex1);
        adjacencyList[v].push_back(sVertex2);

        adjacencyList[u].back().Itr = prev(adjacencyList[v].end());
        adjacencyList[v].back().Itr = prev(adjacencyList[u].end());

        Edge e;
        if(u<v){
            e.U = u;
            e.V = v;
        }
        else {
            e.U = v;
            e.V = u;
        }
        edgeSet.insert(e);
        return;
    }

    void showAdjacencyList()
    {
        for(int i=0; i < adjacencyList.size(); i++)
        {
            list<structureVertex> llist = adjacencyList[i];
//            if(llist.empty())
//                continue;
            for(list<structureVertex>::iterator itr = llist.begin(); itr != llist.end(); itr++)
            {
                structureVertex sVertex = *((*itr).Itr);
                cout << (*itr).Vertex << "->" << sVertex.Vertex << " |";
            }
            cout << "*" <<endl;
        }
    }

    void findDegree()
    {
        structureDegree sDegree;
        // cout << adjacencyList.size() << endl;
        // cout << "started" << endl;
        //Find the degree from the size of the bi-directional adjacency list
        for(unordered_set<int>::iterator itr=vertexSet.begin(); itr != vertexSet.end(); itr++)
        {
            sDegree.Degree = adjacencyList[*itr].size();
            sDegree.Vertex = *itr;
            // cout << "d->"<<sDegree.degree <<"v->"<< sDegree.vertex << endl;
            degreeVector.push_back(sDegree);
        }

        sort(degreeVector.begin(), degreeVector.end(), opDeg);
        // cout << "size" << degreeVector.size() << endl;

        for(int i=0; i<degreeVector.size(); i++)
        {
//            cout << degreeVector[i].Vertex <<" -> " << degreeVector[i].Degree << endl;
            sDegree = degreeVector[i];
            degreeIdxVector[sDegree.Vertex] = i;
        }
//        cout << "---------degreeVector----------"<< endl;
//        print(degreeVector);
//        cout << "---------degreeIdxVector------"<< endl;
//        print(degreeIdxVector);
//        cout << "exited" << endl;
        return;
    }

    void removeEdge(int u, int v)
    {
        if(u > v)
            swap(u,v);
        Edge e;
        e.U = u;
        e.V = v;
        edgeSet.erase(e);
        return;
    }

    void removeVertex(int u)
    {
        int v,idx;
        list<structureVertex>::iterator uitr, vitr;
        Edge e;
        vertexSet.erase(u);

        list<structureVertex>::iterator itr=adjacencyList[u].begin();
        while(itr!=adjacencyList[u].end())
        {
             v = (*itr).Vertex;
             vitr = (*itr).Itr;
             idx = degreeIdxVector[v];
             removeEdge(u, v);

//             cout << "vertex to be modified is" << v << endl;

             edgeSet.erase(e);
//             cout << "edge deleted from the edge set" << endl;
             itr = adjacencyList[u].erase(itr);
             adjacencyList[v].erase(vitr);
        }

    }

    void removeVertexBuss(int u)
    {

        structureDegree sDegree;
        int v,idx;
        list<structureVertex>::iterator uitr, vitr;
        Edge e;
        vertexSet.erase(u);

        idx = degreeIdxVector[u];
        degreeVector[idx].Degree = -1;

        while (idx+1 < degreeVector.size() && degreeVector[idx].Degree < degreeVector[idx+1].Degree)
        {
//                 cout << "gone inside-?" << idx << endl;
            swap(degreeVector[idx], degreeVector[idx+1]);
            swap(degreeIdxVector[degreeVector[idx].Vertex], degreeIdxVector[degreeVector[idx+1].Vertex]);
            idx++;
        }


//        cout << "vertex deleted " << u << "from the vertex set" << endl;

        list<structureVertex>::iterator itr=adjacencyList[u].begin();
        while(itr!=adjacencyList[u].end())
        {
             v = (*itr).Vertex;
             vitr = (*itr).Itr;
             idx = degreeIdxVector[v];

             removeEdge(u,v);
//             cout << "edge deleted from the edge set" << endl;
             itr = adjacencyList[u].erase(itr);
             adjacencyList[v].erase(vitr);


             sDegree = degreeVector[idx];
             degreeVector[idx].Degree = sDegree.Degree - 1;

//             cout << "degree reduction done for vertex" << v << endl;
//             cout << "idx ->" << idx << "size-> " << degreeVector.size() << endl;
             if(degreeVector[idx].Degree == 0)
                 vertexSet.erase(v);

             while (idx+1 < degreeVector.size() && degreeVector[idx].Degree < degreeVector[idx+1].Degree)
             {
//                 cout << "gone inside-?" << idx << endl;
                 swap(degreeVector[idx], degreeVector[idx+1]);
                 swap(degreeIdxVector[degreeVector[idx].Vertex], degreeIdxVector[degreeVector[idx+1].Vertex]);
                 idx++;
             }


//             cout << "---------degreeVector----------"<< endl;
//             print(degreeVector);
//             cout << "---------degreeIdxVector------"<< endl;
//             print(degreeIdxVector);

//             cout << "printing the adjacency list" << endl;
//             showAdjacencyList();
//             cout << "sort the degree finished for vertex " << v << endl;
//             print(degreeVector);
        }

    }



    bool bussReduction(int& k, vector<int>& necessaryVC)
    {
//        cout << "buss reduction" << endl;
        int initialK = k;
        degreeVector.clear();
        findDegree();

        int j = degreeVector.size()-1;

        //remove low degree vertices
        while(j>-1 && degreeVector[j].Degree <= 0)
        {
            vertexSet.erase(degreeVector[j].Vertex);
            j--;
        }

        int i=0;
        int hdVertexCount=0;
        //remove high degree vertices

        while(vertexSet.size() > 0 && degreeVector[0].Degree > k)
        {
            necessaryVC.push_back(degreeVector[i].Vertex);
            removeVertexBuss(degreeVector[i].Vertex);
            k--;
        }

//        cout << "hdcount ->"<<hdVertexCount << " vertex set size ->" << vertexSet.size() << " edge set size" << edgeSet.size() << endl;

        //check the size of necessaryVertex if it is more than k then return false
        if(k < 0)
            return false;
        if( vertexSet.size() > initialK *(initialK +1) || edgeSet.size() > initialK * initialK)
            return false;

        return true;

    }
    bool driverBusReduction(int k, vector<int>& necessaryVC, vector<int>& potentialVC)
    {
        int flag = 2;


        return true;

//        cout << "*** k->"<<k<<"** flag->"<<flag<<"***" <<endl;
//        cout << "printing the necessary vertex" <<endl;
//        print(necessaryVC);

    }

    bool verifierVC(vector<int>& potentialVCVector, unordered_set<Edge, EdgeHashFunction>& edgeSetCheck)
    {
        int v1, v2;
        unordered_set<int> potentialVCSet;

        for(int j=0; j < potentialVCVector.size(); j++)
            potentialVCSet.insert(potentialVCVector[j]);

        unordered_set<Edge, EdgeHashFunction> :: iterator itr;
        Edge e;
        for(itr = edgeSetCheck.begin(); itr != edgeSetCheck.end(); itr++)
        {
           e = *itr;
           v1 = e.U;
           v2 = e.V;
           if(potentialVCSet.find(v1)==potentialVCSet.end()  && potentialVCSet.find(v2)==potentialVCSet.end())
               return false;
        }
        return true;
    }


    bool generateVCPowerSetMethod(vector<int> vertexList, unordered_set<Edge, EdgeHashFunction> curEdgeSet, int k, vector<int>& potentialVC)
    {
        int n = vertexList.size();
        long long int pow_n = (int) pow(2, n);
        int count;
//        cout << "k is" << k <<" "<< pow_n << endl;

//        print(vertexList);
        // Run counter i from 000..0 to 111..1
        for (int i = 1; i < pow_n ; i++)
        {

            potentialVC.clear();
            // consider each element in the set
            count = k;
            for (int j = 0; j < n; j++)
            {
                if ((i & (1 << j)) != 0)
                {
                    potentialVC.push_back(vertexList[j]);
                    count--;
                }
            }
//            cout << "came" << endl;
//            print(potentialVC);
            if(count >=0 && verifierVC(potentialVC, curEdgeSet))
                return true;
        }
        return false;

    }

    //use greedy algorithm for finding maximal matching
    void findMaximalMatching()
    {
        maximalMatching.clear();
        unordered_set<int> maximalvertexSet;
        Edge e;
        for(unordered_set<Edge, EdgeHashFunction>::iterator itr = edgeSet.begin(); itr != edgeSet.end(); itr++)
        {
            e = *itr;
            if(maximalvertexSet.find(e.U) != maximalvertexSet.end() || maximalvertexSet.find(e.V) != maximalvertexSet.end())
                continue;
            maximalvertexSet.insert(e.U);
            maximalvertexSet.insert(e.V);
            maximalMatching.push_back(e);
        }
//        cout << "max matching size -> " << maximalMatching.size() << endl;
//        print(maximalMatching);
        return;
    }


    
    bool findAugmentingPathBPGraph(vector<int>& augmentingPath, vector<int>biPartiteX, vector<int>biPartiteY, vector<Edge> Matching, unordered_map<int, list<int>> adjacencyListBPMap)
    {
        bool augpathFlag = false;

        unordered_map<int, bool> saturatedDic;
        unordered_map<int, bool> markDic;
        unordered_map<int, int> parentDic;
        unordered_map<int, int> matchEdgeDic;

//        cout << "in the augmenting path function" << endl;

        for(int i=0; i<biPartiteX.size(); i++)
        {
            markDic.insert({biPartiteX[i], false});
            saturatedDic.insert({biPartiteX[i], false});
            parentDic.insert({biPartiteX[i], -1});
        }

        for(int i=0; i<biPartiteY.size(); i++)
        {
            markDic.insert({biPartiteY[i], false});
            saturatedDic.insert({biPartiteY[i], false});
            parentDic.insert({biPartiteY[i], -1});
        }

        for(int i=0; i<Matching.size(); i++)
        {
            saturatedDic[Matching[i].U] = true;
            saturatedDic[Matching[i].V] = true;

            matchEdgeDic.insert({Matching[i].U, Matching[i].V});
            matchEdgeDic.insert({Matching[i].V, Matching[i].U});
        }

        for(int i=0; i < biPartiteX.size(); i++)
            if(saturatedDic[biPartiteX[i]] == false)
                augpathFlag = true;

        if(!augpathFlag)
            return augpathFlag;


        queue<int> bfsQueue;
        int v, u, node;
        for(int i=0; i < biPartiteX.size(); i++)
        {
            if(saturatedDic[biPartiteX[i]])
                continue;
            if(!markDic[biPartiteX[i]])
                bfsQueue.push(biPartiteX[i]);
            while(!bfsQueue.empty()){
//                cout << "bfsqueue while loop " <<endl;
              v = bfsQueue.front();
              bfsQueue.pop();
              markDic[v] = true;
//              cout << "in the for loop with vertex -> "<< v << " <-" << endl;
              list<int> neighbourV = adjacencyListBPMap[v];
              for(list<int>::iterator itr = neighbourV.begin(); itr != neighbourV.end(); itr++)
              {
                 u = *itr;
                 if(markDic[u])
                     continue;
                 parentDic[u] = v;
//                 cout << "neighbouring for loop with neighbour vertex -> " << u << " <-" << endl;
                 if(!saturatedDic[u])
                 {
                     node = u;
                     while(node != -1)
                     {
                         augmentingPath.push_back(node);
                         node = parentDic[node];
                     }

                     return true;
                 }
                 else {
                     markDic[u] = true;
                     parentDic[matchEdgeDic[u]] = u;
                     bfsQueue.push(matchEdgeDic[u]);
                 }
              }
            }
        }
        return false;

    }


    
    void findMaximumMatching(vector<int>biPartiteX, vector<int>biPartiteY, unordered_map<int, list<int>>adjacencyListBPMap, vector<Edge>& matching)
    {

        vector<Edge> newMatching;
        vector<int> augmentingPath;
        while(findAugmentingPathBPGraph(augmentingPath, biPartiteX, biPartiteY, matching, adjacencyListBPMap))
        {
            newMatching.clear();
            delta(augmentingPath, matching, newMatching);
            matching = newMatching;
            augmentingPath.clear();
        }

//        cout << "***************the maximum matching for the bipartite graph is shown below**********" << endl;
//        print(matching);
        return;
    }


    void findMinimumVertexCoverFromMaximumMatching(vector<int> bipartiteX, vector<int> bipartiteY, vector<Edge>& matching, vector<int>& minVertexCover, unordered_set<Edge, EdgeHashFunction> bpEdgeSet)
    {
        unordered_map<int, list<int>> adjacencyListBPMap;
        Edge edge;
        unordered_set<int> bipartiteSetX;
        unordered_set<int> bipartiteSetY;
        unordered_set<Edge, EdgeHashFunction> matchingSet;

        for(int i=0; i<matching.size(); i++)
            matchingSet.insert(matching[i]);

        for(int i=0; i<bipartiteX.size(); i++)
            bipartiteSetX.insert(bipartiteX[i]);
        for(int i=0; i<bipartiteY.size(); i++)
            bipartiteSetY.insert(bipartiteY[i]);

        for(int i=0; i<matching.size(); i++)
        {
            edge = matching[i];
            if(bipartiteSetX.find(edge.U) != bipartiteSetX.end())
                adjacencyListBPMap[edge.V].push_back(edge.U);
            else
                adjacencyListBPMap[edge.U].push_back(edge.V);
        }

        for(unordered_set<Edge,EdgeHashFunction>::iterator itr = bpEdgeSet.begin(); itr != bpEdgeSet.end(); itr++)
        {
            edge = *itr;
            if(matchingSet.find(edge) != matchingSet.end())
                continue;
            if(bipartiteSetX.find(edge.U) != bipartiteSetX.end())
                adjacencyListBPMap[edge.U].push_back(edge.V);
            else
                adjacencyListBPMap[edge.V].push_back(edge.U);
        }

//        for(unordered_map<int, list<int>>::iterator itr = adjacencyListBPMap.begin(); itr != adjacencyListBPMap.end(); itr++)
//        {
//            cout << itr->first << " -> ";
//            list<int> llist = itr->second;

//            for(list<int>::iterator itr = llist.begin(); itr != llist.end(); itr++)
//            {
//                int vertex = *itr;
//                cout << vertex  << " |";
//            }
//            cout << "*" <<endl;
//        }


        unordered_map<int, bool> saturatedDic;
        unordered_map<int, bool> markDic;

        for(int i=0; i<bipartiteX.size(); i++)
        {
            saturatedDic.insert({bipartiteX[i], false});
            markDic.insert({bipartiteX[i], false});
        }

        for(int i=0; i<bipartiteY.size(); i++)
        {
            saturatedDic.insert({bipartiteY[i], false});
            markDic.insert({bipartiteY[i], false});
        }

        for(int i=0; i<matching.size(); i++)
        {
            edge = matching[i];
            saturatedDic[edge.U] = true;
            saturatedDic[edge.V] = true;
        }

        queue<int> bfsQueue;
        int u, v;


        for(int i=0; i<bipartiteX.size(); i++)
        {
            if(saturatedDic[bipartiteX[i]])
                continue;
            if(!markDic[bipartiteX[i]])
                bfsQueue.push(bipartiteX[i]);

            while(!bfsQueue.empty())
            {
                v = bfsQueue.front();
                bfsQueue.pop();
                markDic[v] = true;
                list<int> llist = adjacencyListBPMap[v];
                for(list<int>::iterator itr = llist.begin(); itr!=llist.end(); itr++)
                {
                    if(markDic[*itr])
                        continue;
                    bfsQueue.push(*itr);
                }

            }

        }
        vector<int> zVector;
        for(unordered_map<int, bool>::iterator itr = markDic.begin(); itr != markDic.end(); itr++)
        {
            if(itr->second)
            {
//                cout << itr->first << endl;
                zVector.push_back(itr->first);
            }
        }

        sort(bipartiteX.begin(), bipartiteX.end());
        sort(zVector.begin(), zVector.end());
        sort(bipartiteY.begin(), bipartiteY.end());
//        print(bipartiteX);
//        print(zVector);
//        print(bipartiteY);


        int i = 0, j = 0;
        while(i != bipartiteX.size() && j != zVector.size())
        {
            if(bipartiteX[i] == zVector[j])
            {
                i++;
                j++;
            }
            else if(bipartiteX[i]  < zVector[j])
            {
                minVertexCover.push_back(bipartiteX[i]);
                i++;
            }
            else {
                j++;
            }
        }

        if(j == zVector.size()){
            while(i != bipartiteX.size())
            {
                minVertexCover.push_back(bipartiteX[i]);
                i++;
            }
        }

        i = 0;
        j = 0;

        while(i != bipartiteY.size() && j != zVector.size())
        {
            if(bipartiteY[i] == zVector[j])
            {
                minVertexCover.push_back(bipartiteY[i]);
                i++;
                j++;
            }
            else if(bipartiteY[i]  < zVector[j])
                i++;
            else {
                j++;
            }
        }

//        cout << "******************the minimum vertex cover for the given graph is shown below************" << endl;
//        print(minVertexCover);



    }




    bool singleDegreeReductionRule(vector<int>& necessaryCover, int& k)
    {
//        cout << "single degree reduction rule" << endl;
        if(degreeVector.empty())
            return true;
        vector<structureDegree>::iterator itr = degreeVector.end()-1;
        list<structureVertex> llist;
        structureVertex svertex;
        structureDegree sdegree;
        while(itr != degreeVector.begin())
        {
            sdegree = *itr;
            if(sdegree.Degree == 1)
            {
                llist = adjacencyList[sdegree.Vertex];
                svertex = llist.front();
                necessaryCover.push_back(svertex.Vertex);
                removeVertexBuss(sdegree.Vertex);
                removeVertexBuss(svertex.Vertex);
                k--;
            }
            if(k < 0)
                return false;

            else if(sdegree.Degree > 1)
                break;
            itr--;
        }
        return true;
    }

    bool isContained(int u, int v, vector<int>& necessaryCover, int& k)
    {
        vector<int> setA, setB, result;
        list<structureVertex> llist;
        llist = adjacencyList[u];
        for(list<structureVertex>::iterator listitr = llist.begin(); listitr != llist.end(); listitr++)
            setA.push_back(listitr->Vertex);
        llist = adjacencyList[v];
        for(list<structureVertex>::iterator listitr = llist.begin(); listitr != llist.end(); listitr++)
            setB.push_back(listitr->Vertex);
        setIntersection(setA, setB, result);
        if(result.size() == setB.size())
        {
            necessaryCover.push_back(u);
            k--;
            return true;
        }
        return false;
    }

    bool neighbourhoodContainingReductionRule(vector<int>& necessaryCover, int& k)
    {
//        cout << "inside neighbourhood" << endl;
        unordered_set<Edge>::iterator itr;
        Edge edge;

        itr = edgeSet.begin();
        while(itr != edgeSet.end())
        {
            edge = *itr;
            if (adjacencyList[edge.U].size() <= adjacencyList[edge.V].size())
                isContained(edge.U, edge.V, necessaryCover, k);
            else
                isContained(edge.V, edge.U, necessaryCover, k);
            itr++;
        }
        if(k>=0)
            return true;
        else {
            return false;
        }
    }


    bool checkUniqueOptimumHalfIntegralSolution(unordered_map<int, int>& lpDic, int prevlpSum, vector<int>& necessaryVertexCover, int& k)
    {
        //this function has a temporary storage which stores what are the edges that are deleted
        //during a vertex deletion for the optimality checking and restores back after finding the
        //weight of vertex cover in the lp solution.

        list<structureVertex> tempAdjList;
        unordered_map<int, int>::iterator vertexItr;
        int vertex;
        list<structureVertex>::iterator listItr;
        unordered_map<int, int> nlpDic;
        int currlpSum;

        for(vertexItr = lpDic.begin(); vertexItr != lpDic.end(); vertexItr++)
        {
//            degreeVector.clear();
//            findDegree();

            vertex = vertexItr->first;
            tempAdjList =  adjacencyList[vertex];

            removeVertex(vertex);
            currlpSum = findLPSolution(nlpDic);

            vertexSet.insert(vertex);
            for(listItr = tempAdjList.begin(); listItr != tempAdjList.end(); listItr++)
                addEdge(vertex, listItr->Vertex);

            if(currlpSum <= prevlpSum-2)
            {
                removeVertex(vertex);
                lpDic.erase(vertex);
                necessaryVertexCover.push_back(vertex);
                return false;
            }
        }
        return true;
    }

    int findLPSolution(unordered_map<int, int>& lpDic)
    {
        vector<int> biPartiteX, biPartiteY;
        unordered_set<int>::iterator itr;
        int vertex;
        for(itr = vertexSet.begin(); itr != vertexSet.end(); itr++)
        {
            vertex = *itr;
            biPartiteX.push_back(2*vertex-1);
            biPartiteY.push_back(2*vertex);
        }

        unordered_map<int, list<int>> adjacencyListBPMap;
        unordered_set<Edge, EdgeHashFunction> bpEdgeSet;
        Edge e1, e2;
        unordered_set<Edge>::iterator edgeItr;
        Edge edge;
        for(edgeItr = edgeSet.begin(); edgeItr != edgeSet.end(); edgeItr++)
        {
            edge = *edgeItr;
            if(edge.U < edge.V)
            {
                e1.U = 2*edge.U-1;
                e1.V = 2*edge.V;
                e2.U = 2*edge.U;
                e2.V = 2*edge.V-1;
            }
            else
            {
                e1.U = 2*edge.V-1;
                e1.V = 2*edge.U;
                e2.U = 2*edge.V;
                e2.V = 2*edge.U-1;
            }
            bpEdgeSet.insert(e1);
            bpEdgeSet.insert(e2);
            adjacencyListBPMap[2*edge.U-1].push_back(2*edge.V);
            adjacencyListBPMap[2*edge.V].push_back(2*edge.U-1);
            adjacencyListBPMap[2*edge.V-1].push_back(2*edge.U);
            adjacencyListBPMap[2*edge.U].push_back(2*edge.V-1);
        }

//        cout <<"new edge set size is->"<< bpEdgeSet.size() << endl;

        vector<Edge> matching;
        vector<int> minVertexCover;

        findMaximumMatching(biPartiteX, biPartiteY, adjacencyListBPMap, matching);
        findMinimumVertexCoverFromMaximumMatching(biPartiteX, biPartiteY, matching, minVertexCover, bpEdgeSet);

        unordered_set<int> vertexCoverSet;
        for(int i=0; i<minVertexCover.size(); i++)
            vertexCoverSet.insert(minVertexCover[i]);

        int lpSum=0;

        for(itr = vertexSet.begin(); itr!=vertexSet.end(); itr++)
        {
            vertex = *itr;
            if(vertexCoverSet.find(2*vertex-1) != vertexCoverSet.end() && vertexCoverSet.find(2*vertex) != vertexCoverSet.end())
            {
                lpDic[vertex] = 2;
                lpSum += 2;
            }
            else if(vertexCoverSet.find(2*vertex-1) != vertexCoverSet.end() || vertexCoverSet.find(2*vertex) != vertexCoverSet.end())
            {
                lpDic[vertex] = 1;
                lpSum += 1;
            }
            else
                lpDic[vertex] = 0;
        }
        return lpSum;
    }


    bool randomizedBranchingAlgorithm(vector<int>& potentialVC, int k)
    {
//        cout << "inside random branching" << endl;
//        cout << "inside branching" << endl;
        if(k < 0)
            return false;

        if(k == 0 && !edgeSet.empty())
            return false;

        if(edgeSet.empty())
            return true;

        if(vertexSet.empty())
            return false;


        int randomBit = rand()%2;

        unordered_set<int>::iterator vtxItr;
        vtxItr = vertexSet.begin();
        int vertex = *vtxItr;


        vector<int> neighbor;
        vector<Edge> removedEdges;
        vector<Edge> removedNeighbourEdges;
        Edge edge;

        list<structureVertex> nlist = adjacencyList[vertex];


        for(list<structureVertex>::iterator itr=nlist.begin(); itr != nlist.end(); itr++)
            neighbor.push_back(itr->Vertex);

        edge.U = vertex;
        for(int i=0; i < neighbor.size(); i++)
        {
            edge.V = neighbor[i];
            removedEdges.push_back(edge);
        }

        removeVertex(vertex);

        int nvertex;

        for(int i = 0; i < neighbor.size(); i++){
            nvertex = neighbor[i];
            nlist = adjacencyList[nvertex];
            edge.U = nvertex;
            for(list<structureVertex>::iterator itr=nlist.begin(); itr != nlist.end(); itr++)
            {
                edge.V = itr->Vertex;
                removedNeighbourEdges.push_back(edge);
            }
            removeVertex(nvertex);

        }

        if(randomBit){

            for(int i=0; i<neighbor.size(); i++)
            {
                potentialVC.push_back(neighbor[i]);
            }

            if(randomizedBranchingAlgorithm(potentialVC, k-neighbor.size()))
               return true;

            for(int i=0; i<neighbor.size(); i++)
            {
                vertexSet.insert(neighbor[i]);
                potentialVC.pop_back();
            }

            for(int i=0; i<removedNeighbourEdges.size(); i++)
                addEdge(removedNeighbourEdges[i].U, removedNeighbourEdges[i].V);

            potentialVC.push_back(vertex);

            if(randomizedBranchingAlgorithm(potentialVC, k-1))
                return true;

            potentialVC.pop_back();
            vertexSet.insert(vertex);

            for(int i=0; i<removedEdges.size(); i++)
                addEdge(removedEdges[i].U, removedEdges[i].V);

            return false;
        }
        else {

            for(int i=0; i<neighbor.size(); i++)
            {
                vertexSet.insert(neighbor[i]);
            }

            for(int i=0; i<removedNeighbourEdges.size(); i++)
                addEdge(removedNeighbourEdges[i].U, removedNeighbourEdges[i].V);

            potentialVC.push_back(vertex);

            if(randomizedBranchingAlgorithm(potentialVC, k-1))
                return true;

            potentialVC.pop_back();

            for(int i=0; i<neighbor.size(); i++)
            {
                potentialVC.push_back(neighbor[i]);
                removeVertex(neighbor[i]);
            }

            if(randomizedBranchingAlgorithm(potentialVC, k-neighbor.size()))
               return true;

            vertexSet.insert(vertex);

            for(int i=0; i<neighbor.size(); i++)
            {
                vertexSet.insert(neighbor[i]);
                potentialVC.pop_back();
            }

            for(int i=0; i<removedEdges.size(); i++)
                addEdge(removedEdges[i].U, removedEdges[i].V);

            for(int i=0; i<removedNeighbourEdges.size(); i++)
                addEdge(removedNeighbourEdges[i].U, removedNeighbourEdges[i].V);

            return false;

        }



    }


    bool branchingAlgorithm(vector<int>& potentialVC, int k)
    {
//        cout << "inside branching" << endl;
        if(k < 0)
            return false;

        if(k == 0 && !edgeSet.empty())
            return false;

        if(edgeSet.empty())
            return true;

        if(vertexSet.empty())
            return false;


        unordered_set<int>::iterator vtxItr;
        vtxItr = vertexSet.begin();
        int vertex = *vtxItr;

        potentialVC.push_back(vertex);
        vector<int> neighbor;
        vector<Edge> removedEdges;
        Edge edge;

        list<structureVertex> nlist = adjacencyList[vertex];


        for(list<structureVertex>::iterator itr=nlist.begin(); itr != nlist.end(); itr++)
            neighbor.push_back(itr->Vertex);

        edge.U = vertex;
        for(list<structureVertex>::iterator itr=nlist.begin(); itr != nlist.end(); itr++)
        {
            edge.V = itr->Vertex;
            removedEdges.push_back(edge);
        }
        removeVertex(vertex);

        if(branchingAlgorithm(potentialVC, k-1))
            return true;

        potentialVC.pop_back();

        int nvertex;
        for(int i = 0; i < neighbor.size(); i++){
            nvertex = neighbor[i];
            nlist = adjacencyList[nvertex];
            edge.U = nvertex;
            for(list<structureVertex>::iterator itr=nlist.begin(); itr != nlist.end(); itr++)
            {
                edge.V = itr->Vertex;
                removedEdges.push_back(edge);
            }
            potentialVC.push_back(nvertex);
            removeVertex(neighbor[i]);
        }


        if(branchingAlgorithm(potentialVC, k-neighbor.size()))
           return true;


        vertexSet.insert(vertex);

        for(int i=0; i<neighbor.size(); i++)
        {
            vertexSet.insert(neighbor[i]);
            potentialVC.pop_back();
        }

        for(int i=0; i<removedEdges.size(); i++)
            addEdge(removedEdges[i].U, removedEdges[i].V);

        return false;
    }



    bool driverlp(vector<int>& necessaryVertexCover, int& k, vector<int>& potentialVertexCover, int randmode)
    {
//        cout << "inside driver lp" << endl;
        unordered_map<int, int> lpDic;
        unordered_set<int>::iterator itr;

        for(itr = vertexSet.begin(); itr!=vertexSet.end(); itr++)
            lpDic[*itr] = -1;

        bool flag = false;
        int lpSum;
        int oldSize;
        int initialK = k;



//        cout << "came in driver lp" << endl;


        while(!lpDic.empty() && k>0)
        {
            oldSize = lpDic.size();
            lpSum = findLPSolution(lpDic);
            if(lpSum > 2*k)
                return false;

            unordered_map<int, int>::iterator lpitr;
            int key, value;

            lpitr = lpDic.begin();
            while( lpitr != lpDic.end())
            {
                key = lpitr->first;
                value = lpitr->second;
                if(value == 0)
                {
                    removeVertex(key);
                    lpitr = lpDic.erase(lpitr);
                }
                else if(value == 2)
                {
                    removeVertex(key);
                    lpitr = lpDic.erase(lpitr);
                    necessaryVertexCover.push_back(key);
                    k--;
                }
                else {
                    lpitr++;
                }

            }

            if( oldSize == lpDic.size() && checkUniqueOptimumHalfIntegralSolution(lpDic, lpSum, necessaryVertexCover, k) )
            {
                flag = true;
                break;
            }

        }

        if(lpDic.size() > 2*initialK)
            return false;


        unordered_set<int> vertexSetR;
        unordered_map<int, bool> visit;
        vertexSetR = vertexSet;

//        cout << "after lp "<<"vertex size-> " << vertexSet.size() << " edge size -> " << edgeSet.size() << endl;

//        cout << "going to recursion" << endl;
        if(flag)
        {
            if(!randmode)
                return branchingAlgorithm(potentialVertexCover, k);
            else {
                return randomizedBranchingAlgorithm(potentialVertexCover, k);
            }
        }
        else if(lpDic.empty())
            return true;
        else
            return false;


    }







    bool FindCrown(int k, vector<int> head)
    {
//        cout << "inside crown" << endl;
        findMaximalMatching();
        if(maximalMatching.size() > k)
            return false;

//        cout << "maximal match passed" << endl;
//        print(maximalMatching);

        vector<int> Vm;
        unordered_set<int> VmSet;
        vector<int> I;
        unordered_set<int> ISet;

        for(int i=0; i < maximalMatching.size(); i++)
        {
            Vm.push_back(maximalMatching[i].U);
            Vm.push_back(maximalMatching[i].V);
        }

//        Vm.clear();
//        Vm = {5, 7, 6};
        for(int i=0; i<Vm.size(); i++)
        {
            VmSet.insert(Vm[i]);
        }

//        cout << "vertex in maximal matching" << endl;
//        print(Vm);



        vector<int> vertexVector;

        for(unordered_set<int>::iterator itr = vertexSet.begin(); itr != vertexSet.end(); itr++)
            vertexVector.push_back(*itr);

        setMinus(vertexVector, Vm, I);
//        I.clear();
//        I = {1, 2, 3, 4};
//        cout << "independent set" << endl;
//        print(I);


        for(int i=0; i<I.size(); i++)
            ISet.insert(I[i]);


        unordered_map<int, list<int>> adjacencyListBPMap;
        Edge edge;
        unordered_set<Edge, EdgeHashFunction> bpEdgeSet;

        for(unordered_set<Edge>::iterator itr = edgeSet.begin(); itr != edgeSet.end(); itr++)
        {
            edge = *itr;
            if( (VmSet.find(edge.U) != VmSet.end() && ISet.find(edge.V) != ISet.end() ) || (VmSet.find(edge.V) != VmSet.end() && ISet.find(edge.U) != ISet.end() ) )
            {
//                cout << edge.U << "-" << edge.V << endl;
                bpEdgeSet.insert(edge);
                adjacencyListBPMap[edge.U].push_back(edge.V);
                adjacencyListBPMap[edge.V].push_back(edge.U);
            }
        }

        vector<Edge> matching;
        vector<int> vertexCover;

        findMaximumMatching(Vm, I, adjacencyListBPMap, matching);
//        cout << "maximum matching is done and the result is given" << endl;
//        print(matching);


        if(matching.size() > k)
            return false;

        findMinimumVertexCoverFromMaximumMatching(Vm, I, matching, vertexCover, bpEdgeSet);
//        cout << "the minimum vertex cover is given below" << endl;
//        print(vertexCover);



        setIntersection(vertexCover, Vm, head);

//        cout << "head part of the graph" << endl;
//        print(head);

//        unordered_set<int> HSet;

//        for(int i=0; i<head.size(); i++)
//            HSet.insert(head[i]);


//        for(int i=0; i < matching.size(); i++)
//        {
//            edge = matching[i];
//            if( HSet.find(edge.U) != HSet.end() && ISet.find(edge.V) != ISet.end() )
//                crown.push_back(edge.V);
//            else if( HSet.find(edge.V) != HSet.end() && ISet.find(edge.U) != ISet.end() )
//                crown.push_back(edge.U);
//        }

//        cout << "crown is given below" << endl;
//        print(crown);
//        cout << "end of one iteration" << endl;

    }


    bool crownReduction(int& k, vector<int>& necessaryVC)
    {
//        cout << "inside crown reduction" << endl;
        vector<int> head;
        if(vertexSet.size() < 3*k)
            return true;
        if(FindCrown(k, head))
        {
            for(int i=0; i<head.size(); i++)
            {
                removeVertex(head[i]);
                necessaryVC.push_back(head[i]);
                k--;
            }
            if(k<0)
                return false;
        }
        return true;
    }


    bool findingVCMethod1(int k, vector<int>& potentialVC)
    {
        vector<int> vertexList;

        for(unordered_set<int>::iterator itr = vertexSet.begin(); itr != vertexSet.end(); itr++)
            vertexList.push_back(*itr);
        if ( generateVCPowerSetMethod(vertexList, edgeSet, k, potentialVC))
        {
//            cout << "with k value " << k << "vertex cover is" << endl;
//            print(potentialVC);
            return true;
        }
        else {
//            cout << "no vertex cover of size" << k << "exist" <<endl;
            return false;
        }
    }

    bool findingVCMethod2(int k, vector<int>& potentialVC)
    {
        vector<int> vertexVector;
        for(unordered_set<int>::iterator itr = vertexSet.begin(); itr != vertexSet.end(); itr++){
            vertexVector.push_back(*itr);
        }
        vector<vector<int>> ans;

        ans = makeCombi(k, vertexVector);

        for(int i=0; i<ans.size(); i++)
        {
            if(verifierVC(ans[i], edgeSet)){
                potentialVC = ans[i];
//                cout << "vertex cover found of size" << k << endl;
//                print(potentialVC);
                return true;
            }
        }
//        cout << "no vertex cover of size" << k << "found" << endl;
        return false;

    }

    bool findingVCMethod3(int k, vector<int> necessaryVC, vector<int>& potentialVC)
    {
        unordered_set<int> vertexSetR;
        vertexSetR = vertexSet;
        unordered_map<int, bool> visit;
        if(branchingAlgorithm(potentialVC, k))
        {
//            print(potentialVC);
            return true;
        }
        else {
            return false;
        }

    }

    bool findingVCMethod4(int k, vector<int>& necessaryVC, vector<int>& potentialVC)
    {
        //buss reduction
//        cout << "initial "<<"vertex size-> " << vertexSet.size() << " edge size -> " << edgeSet.size() << endl;


        int oldk = INT_MAX;
        while(oldk > k)
        {
//            cout <<"**"<< k << endl;
            oldk = k;
            if(!bussReduction(k, necessaryVC))
            {
//                cout << "no vertex cover of that size found"<< endl;
                return false;
            }
        }

//        cout << "after buss "<<"vertex size-> " << vertexSet.size() << " edge size -> " << edgeSet.size() << endl;


        if(edgeSet.empty())
            return true;

        if(!singleDegreeReductionRule(necessaryVC, k)){
//            cout << "no vertex cover of that size found";
            return false;
        }

//        cout << "after single "<<"vertex size-> " << vertexSet.size() << " edge size -> " << edgeSet.size() << endl;


        if(edgeSet.empty())
            return true;

        if(!neighbourhoodContainingReductionRule(necessaryVC, k))
        {
//                cout << "no vertex cover of that size found";
                return false;
        }

//        cout << "after containing "<<"vertex size-> " << vertexSet.size() << " edge size -> " << edgeSet.size() << endl;


        if(edgeSet.empty())
            return true;


        oldk = k+1;
        while(oldk > k)
        {
            oldk = k;
            if(!crownReduction(k, necessaryVC))
            {
//                cout << "no vertex cover of that size found" << endl;
                return false;
            }
            if(!bussReduction(k, necessaryVC))
                return false;
        }

//        cout << "after crown "<<"vertex size-> " << vertexSet.size() << " edge size -> " << edgeSet.size() << endl;


        if(edgeSet.empty())
            return true;

        if(driverlp(necessaryVC, k, potentialVC, 0))
          return true;
        else
            return false;

    }

    bool findingVCMethod5(int k, vector<int>& necessaryVC, vector<int>& potentialVC)
    {
        int oldk = INT_MAX;
        while(oldk > k)
        {
//            cout <<"**"<< k << endl;
            oldk = k;
            if(!bussReduction(k, necessaryVC))
            {
//                cout << "no vertex cover of that size found"<< endl;
                return false;
            }
        }
        if(!singleDegreeReductionRule(necessaryVC, k)){
//            cout << "no vertex cover of that size found";
            return false;
        }

        oldk = k+1;
        while(oldk > k)
        {
            oldk = k;
            if(!crownReduction(k, necessaryVC))
            {
//                cout << "no vertex cover of that size found" << endl;
                return false;
            }
            if(!bussReduction(k, necessaryVC))
                return false;
        }

        if(driverlp(necessaryVC, k, potentialVC, 1))
          return true;
        else
            return false;

    }

    bool findingVCMethod6(int k, vector<int> necessaryVC, vector<int>& potentialVC)
    {
        unordered_set<int> vertexSetR;
        vertexSetR = vertexSet;
        unordered_map<int, bool> visit;
        if(randomizedBranchingAlgorithm(potentialVC, k))
        {
//            print(potentialVC);
            return true;
        }
        else {
            return false;
        }

    }


};

void findVertexCover(int nodeNum, int edgeNum, vector<Edge>& edgeset)
{
    Graph graph(nodeNum, edgeNum);
    for(int i=0; i<edgeset.size(); i++){
        graph.addEdge(edgeset[i].U, edgeset[i].V);
    }

    graph.findMaximalMatching();
    int l = graph.maximalMatching.size();
    int h = 2*l;
    int i=0;
    int mid = 0;
    vector<int> necessaryVC;
    vector<int> potentialVC;
//    cout << "it is sadn l "<< l << "and h is " << h << endl;
    while(l < h){
        i++;
        Graph graph(nodeNum, edgeNum);
        for(int i=0; i<edgeset.size(); i++){
            graph.addEdge(edgeset[i].U, edgeset[i].V);
        }
//            cout << "iteration number" << i << endl;
        mid = (l+h)/2;
//        cout << "for the value l "<< l << "and h " << h <<"and the mid value "<<  mid << " the vertex cover is"  << endl;

        necessaryVC.clear();
        potentialVC.clear();


        //this part is to test the buss reduction
        if (graph.findingVCMethod1(mid, potentialVC)){
            h = mid;
        }

        else {
            l = mid+1;
        }

//        print(necessaryVC);
//        print(potentialVC);


    }
    necessaryVC.clear();
    potentialVC.clear();
    mid = (l+h)/2;
//    cout << "for the value l "<< l << "and h " << h <<"and the mid value "<<  mid << "the vertex cover is"  << endl;
    //this part is to test the buss reduction


    if (graph.findingVCMethod1(mid, potentialVC) ){
//            print(necessaryVC);
//        cout << "finally returning the value" << endl;
//        cout << potentialVC.size() << endl;
//        print(necessaryVC);
//        print(potentialVC);
    }
    else {
        cout << "wrong function" << endl;
    }
}



void fileReadingIBM(string path, vector<int>& nodeVector, vector<int>& edgeVector, vector<vector<Edge>>& edgeSetVector)
{
    string mystr;
    ifstream inFile;
    Edge edge;
    int nodeNum = -1;
    int edgeNum = -1;
    vector<Edge> edgeSet;
    int flag = -1;
    int x;

    inFile.open(path);

    if(!inFile)
        cout << "could not open" << endl;
    else
        cout << "opened successfully" << endl;
    int i=0, data, j=0;
    string str;

    while (inFile) {
        // Read a Line from File
        getline(inFile, mystr);
        stringstream ss(mystr);
        i++;
        j = 0;
        while(ss >> str)
        {
            j++;
            if(j == 1){

                if(str == "t")
                {
//                    cout <<"node -> " << nodeNum << " edge-> " << edgeNum << " and the edgeSet is" <<endl;
//                    print(edgeSet);
//                    cin >> x;
                    nodeVector.push_back(nodeNum);
                    edgeVector.push_back(edgeNum);
                    edgeSetVector.push_back(edgeSet);
                    nodeNum = 0;
                    edgeNum = 0;
                    edgeSet.clear();
                    flag = 0;
                }
                else if(str == "v"){
                    flag = 1;
                    nodeNum++;
                }

                else if(str == "e")
                {
                    flag = 2;
                    edgeNum++;
                }
            }

            if(flag == 2){

                if(j == 2)
                    edge.U = stoi(str)+1;
                else if(j == 3)
                    edge.V = stoi(str)+1;
            }
        }
        if(flag == 2)
            edgeSet.push_back(edge);
    }

}





void fileReadingPace(string path, int& nodeNum, int& edgeNum, vector<Edge>& edgeSet)
{
    string mystr;
    ifstream inFile;
    Edge edge;

    inFile.open(path);

    if(!inFile)
        cout << "could not open" << endl;
    else
        cout << "opened successfully" << endl;

    int i=0, data, j=0;
    string str;

    while (inFile) {
        // Read a Line from File
        getline(inFile, mystr);
        stringstream ss(mystr);
        i++;
        j = 0;
        if(i==1){
            while(ss >> str){
                j++;
//                cout << str << " ";
                if(j == 3)
                    nodeNum = stoi(str);
                else if(j == 4)
                    edgeNum = stoi(str);
            }
            continue;
        }
        while(ss >> data)
        {
//            if(i==2)
//                cout << data << endl;
            j++;
            if(j==1)
                edge.U = data;
            else if(j==2)
                edge.V = data;
        }

        edgeSet.push_back(edge);

        // Print line in Console
//         cout << mystr << endl;
    }
    cout <<"---------------------"<< endl << endl << endl;
    cout << "node-> " << nodeNum << " edge-> " << edgeNum << endl;
    for(int i=0; i<edgeSet.size(); i++)
    {
        cout << edgeSet[i].U <<"-----"<< edgeSet[i].V << endl;
    }
    cout << "checked the output" << endl;

}


/*








    bool findVertexCover(reductionType method)
    {
        switch (method)
        {
            case bussReductionMethod:
                bussReduction();

 break;
            case crownReductionMethod:
                crownReduction();
                break;
        }
        return true;
    }


*/




int main()
{
    cout << "Hello World!" << endl;




//*******************this part is for running buss reduction by reading .hgr file***********************


    int nodeNum, edgeNum;
    vector<Edge> edgeSet;
    vector<int> nodeVector;
    vector<int> edgeVector;

    vector<vector<Edge>> edgeSetVector;

    fileReadingIBM("C:\\Users\\Suman\\Documents\\hello_world\\newtest.data", nodeVector, edgeVector, edgeSetVector);
//    cout << "came out of function" << endl;
    int x;
    ofstream myfile ("result_branching_only.txt");
    if (myfile.is_open())
        cout << "opened successfully" << endl;
    else {
        cout << "failed to open" << endl;
    }

    for(int i=1; i<1500; i++)
    {
        nodeNum = nodeVector[i];
        edgeNum = edgeVector[i];
        edgeSet = edgeSetVector[i];
        cout << "calling the funciton" << endl;

//         Get starting timepoint
        auto start = high_resolution_clock::now();
        findVertexCover(nodeNum, edgeNum, edgeSet);
        // Get ending timepoint
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

//          myfile << "This is a line.\n";
//          myfile << "This is another line.\n";

        myfile << i <<" " << duration.count() << endl;
//        cout << "waiting for the user input the last file was numbered "<< i << endl;
//        cin >> x;

    }
    cout << "waiting for the user input the last" << endl;
    cin >> x;
    myfile.close();





//*******************************************************************************************************
//    Graph graph1(nodeNum, edgeNum);
//    for(int i=0; i<edgeSet.size(); i++){
//        graph1.addEdge(edgeSet[i].U, edgeSet[i].V);
//    }
//    vector<int> necessaryvc;
//    vector<int> potentialvc;
//    int k = 100;
//    graph1.driverlp(necessaryvc, k, potentialvc);


//****************this part is to check the symmetric difference fucntion**********************


//    vector<int> augpath = {1, 2, 3, 4, 5, 6};

//    Edge edge;
//    vector<Edge> matching;
//    vector<Edge> newMatching;

//    edge.U = 1, edge.V = 2;
//    matching.push_back(edge);
//    edge.U = 3, edge.V = 4;
//    matching.push_back(edge);

//    delta(augpath, matching, newMatching);

//    cout << " --------------new maximum matching is -------------- " << endl;

//    for(int i=0; i<newMatching.size(); i++)
//        cout << "U -> " << newMatching[i].U << "V -> " << newMatching[i].V << endl;

//**********************************************************************
//        nodeNum = 25;
//        edgeNum = 24;
//        Edge edge;
//        edge.U = 1, edge.V = 2;
//        edgeSet.push_back(edge);
//        edge.U = 1, edge.V = 3;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 4;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 6;
//        edgeSet.push_back(edge);

//        edge.U = 3, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 9;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 10;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 11;
//        edgeSet.push_back(edge);


//        edge.U = 6, edge.V = 12;
//        edgeSet.push_back(edge);
//        edge.U = 6, edge.V = 13;
//        edgeSet.push_back(edge);
//        edge.U = 7, edge.V = 14;
//        edgeSet.push_back(edge);
//        edge.U = 7, edge.V = 15;
//        edgeSet.push_back(edge);
//        edge.U = 8, edge.V = 16;
//        edgeSet.push_back(edge);

//        edge.U = 8, edge.V = 17;
//        edgeSet.push_back(edge);
//        edge.U = 9, edge.V = 18;
//        edgeSet.push_back(edge);
//        edge.U = 9, edge.V = 19;
//        edgeSet.push_back(edge);
//        edge.U = 10, edge.V = 20;
//        edgeSet.push_back(edge);
//        edge.U = 10, edge.V = 21;
//        edgeSet.push_back(edge);

//        edge.U = 11, edge.V = 22;
//        edgeSet.push_back(edge);
//        edge.U = 11, edge.V = 23;
//        edgeSet.push_back(edge);
//        edge.U = 12, edge.V = 24;
//        edgeSet.push_back(edge);
//        edge.U = 12, edge.V = 25;
//        edgeSet.push_back(edge);












//***********************************************************************************************************
//        nodeNum = 6;
//        edgeNum = 5;
//        Edge edge;
//        edge.U = 1, edge.V = 4;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 4;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 4;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 6;
//        edgeSet.push_back(edge);


//****************************************************************************************************
//        nodeNum = 11;
//        edgeNum = 19;
//        Edge edge;
//        edge.U = 1, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 9;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 1, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 6, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 7, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 7, edge.V = 11;
//        edgeSet.push_back(edge);
//        edge.U = 1, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 8, edge.V = 10;
//        edgeSet.push_back(edge);
//        edge.U = 9, edge.V = 10;
//        edgeSet.push_back(edge);
//        edge.U = 10, edge.V = 11;
//        edgeSet.push_back(edge);

//    Graph graph1(11, 19);
//    graph1.addEdge(1,5);
//    graph1.addEdge(5,6);
//    graph1.addEdge(5,8);
//    graph1.addEdge(5,9);
//    graph1.addEdge(2,5);
//    graph1.addEdge(4,5);
//    graph1.addEdge(1,7);
//    graph1.addEdge(3,7);
//    graph1.addEdge(4,7);
//    graph1.addEdge(6,7);
//    graph1.addEdge(7,8);
//    graph1.addEdge(7,11);
//    graph1.addEdge(1,6);
//    graph1.addEdge(2,6);
//    graph1.addEdge(3,6);
//    graph1.addEdge(4,6);
//    graph1.addEdge(8,10);
//    graph1.addEdge(9,10);
//    graph1.addEdge(10,11);


//**************************unique half integral check*************
//    Graph graph1(3, 3);
//    graph1.addEdge(1, 2);
//    graph1.addEdge(2, 3);
//    graph1.addEdge(1, 3);

//    vector<int> necessaryvc;
//    vector<int> potentialvc;
//    int k = 4;

//    graph1.driverlp(necessaryvc, k, potentialvc);

//    cout << k << endl;
//    print(necessaryvc);
//    print(potentialvc);


//    cout << "sum--->" << sum << endl;


//    vector<int> necessaryvc;
//    vector<int> potentialvc;
//    int k = 2;

//    Graph graph1(9, 9);
//    graph1.addEdge(1, 2);
//    graph1.addEdge(2, 3);
//    graph1.addEdge(1, 4);
//    graph1.addEdge(4, 5);
//    graph1.addEdge(2, 5);
//    graph1.addEdge(3, 6);
//    graph1.addEdge(3, 7);
//    graph1.addEdge(3, 8);
//    graph1.addEdge(3, 9);

//    nodeNum = 9;
//    edgeNum = 9;
//    Edge edge;
//    edge.U = 1, edge.V = 2;
//    edgeSet.push_back(edge);
//    edge.U = 2, edge.V = 3;
//    edgeSet.push_back(edge);
//    edge.U = 1, edge.V = 4;
//    edgeSet.push_back(edge);
//    edge.U = 4, edge.V = 5;
//    edgeSet.push_back(edge);
//    edge.U = 2, edge.V = 5;
//    edgeSet.push_back(edge);
//    edge.U = 3, edge.V = 6;
//    edgeSet.push_back(edge);
//    edge.U = 3, edge.V = 7;
//    edgeSet.push_back(edge);
//    edge.U = 3, edge.V = 8;
//    edgeSet.push_back(edge);
//    edge.U = 3, edge.V = 9;
//    edgeSet.push_back(edge);

//    graph1.driverlp(necessaryvc, k, potentialvc);




//    cout << k << endl;
//    print(necessaryvc);
//    print(potentialvc);



//************************************************************************


//        Graph graph1(6, 6);
//        graph1.addEdge(1,4);
//        graph1.addEdge(1,5);
//        graph1.addEdge(1,6);
//        graph1.addEdge(2,5);
//        graph1.addEdge(2,6);
//        graph1.addEdge(3,5);


//        vector<int> biPartiteX = {1, 2, 3};
//        vector<int> biPartiteY = {4, 5, 6};

//***************************************************************************
//        Graph graph1(9, 8);
//        graph1.addEdge(1,6);
//        graph1.addEdge(2,6);
//        graph1.addEdge(2,7);
//        graph1.addEdge(3,8);
//        graph1.addEdge(3,9);
//        graph1.addEdge(4,7);
//        graph1.addEdge(5,6);
//        graph1.addEdge(5,9);

//        nodeNum = 9;
//        edgeNum = 8;
//        Edge edge;
//        edge.U = 1, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 9;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 9;
//        edgeSet.push_back(edge);





//        vector<int> biPartiteX = {1, 2, 3, 4, 5};
//        vector<int> biPartiteY = {6, 7, 8, 9};
//****************************************************************************

//        Graph graph1(8, 7);
//        graph1.addEdge(1,5);
//        graph1.addEdge(1,6);
//        graph1.addEdge(2,7);
//        graph1.addEdge(2,8);
//        graph1.addEdge(3,5);
//        graph1.addEdge(3,8);
//        graph1.addEdge(4,6);

//        nodeNum = 8;
//        edgeNum = 7;
//        Edge edge;
//        edge.U = 1, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 1, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 5;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 6;
//        edgeSet.push_back(edge);


//        vector<int> biPartiteX = {1, 2, 3, 4};
//        vector<int> biPartiteY = {5, 6, 7, 8};
//********************************************************************************

//        Graph graph1(10, 10);
//        graph1.addEdge(1,6);
//        graph1.addEdge(1,10);
//        graph1.addEdge(2,7);
//        graph1.addEdge(2,8);
//        graph1.addEdge(3,6);
//        graph1.addEdge(3,9);
//        graph1.addEdge(4,8);
//        graph1.addEdge(4,9);
//        graph1.addEdge(5,6);
//        graph1.addEdge(5,9);

//        nodeNum = 10;
//        edgeNum = 10;
//        Edge edge;
//        edge.U = 1, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 1, edge.V = 10;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 7;
//        edgeSet.push_back(edge);
//        edge.U = 2, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 3, edge.V = 9;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 8;
//        edgeSet.push_back(edge);
//        edge.U = 4, edge.V = 9;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 6;
//        edgeSet.push_back(edge);
//        edge.U = 5, edge.V = 9;
//        edgeSet.push_back(edge);


//        vector<int> biPartiteX = {1, 2, 3, 4, 5};
//        vector<int> biPartiteY = {6, 7, 8, 9, 10};

//************************************************************************************



//    graph1.findMaximumMatching(biPartiteX, biPartiteY);
//    graph1.findMinimumVertexCoverFromMaximumMatching(biPartiteX, biPartiteY);

//    graph1.addEdge(1,5);

//    graph1.showAdjacencyList();


//    graph1.findDegree();
//    graph1.removeVertexBuss(1);


//    cout << "calling the funciton" << endl;

    // Get starting timepoint
//    auto start = high_resolution_clock::now();
//    findVertexCover(nodeNum, edgeNum, edgeSet);
//    // Get ending timepoint
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);

//    cout << "Time taken by function: "
//         << duration.count() << " microseconds" << endl;





    /*this part is to test the combinatorial generation part
    int k = 3;
    vector<int> data = {1, 2, 4, 5, 6};
    vector<vector<int>> ans = makeCombi(k, data);
    */

    return 0;
}
