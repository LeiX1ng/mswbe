#include "naive.h"
biclique::~biclique()
{
    vertices[0].clear();
    vertices[1].clear();
    left_size=0;
    right_size=0;
    wbc.clear();
    leftMapindex.clear();
    rightMapindex.clear();
    indexMapleft.clear();
    indexMapright.clear();
}
void biclique::get_wbc(biGraph& bg){
    for (auto u : this->vertices[0])
    {
        for (auto v : this->vertices[1])
        {
            float w = bg.edge_weight[make_pair(u, v)];
            weight_biclique temp;
            temp.u=u;
            temp.v=v;
            temp.weight=w;
            this->wbc.push_back(temp);
        }
    }
}
void biclique::print_biclique(){
    cout<<"left size: "<<this->left_size<<endl;
    cout<<"right size: "<<this->right_size<<endl;
    cout<<"left vertices: ";
    for(int i=0;i<this->left_size;i++){
        cout<<this->vertices[0][i]<<" ";
    }
    cout<<endl;
    cout<<"right vertices: ";
    for(int i=0;i<this->right_size;i++){
        cout<<this->vertices[1][i]<<" ";
    }
    cout<<endl;
}
//renumber the vertices in the biclique
void biclique::sorted_wbc(){
    unsigned int leftCount = 0;
    for (const auto &vertex: this->vertices[0]) {
        // 处理左侧顶点
        leftMapindex[vertex] = leftCount;
        indexMapleft[leftCount] = vertex;
        leftCount++;
    }
    unsigned int rightCount = leftCount;
    for (const auto &vertex: this->vertices[1]) {
        // 处理右侧顶点
        rightMapindex[vertex] = rightCount;
        indexMapright[rightCount] = vertex;
        rightCount++;
    }
    for(auto& edge: this->wbc){
        edge.u = leftMapindex[edge.u];
        edge.v = rightMapindex[edge.v];
    }
}
void biclique::print_wbc(){
    cout<<"left size: "<<this->left_size<<endl;
    cout<<"right size: "<<this->right_size<<endl;
    cout<<"weight edges: "<<endl;
    for(auto edge:this->wbc){
        cout<<edge.u<<" "<<edge.v<<" "<<edge.weight<<endl;
    }
}
void printResult(const std::vector<std::pair<std::vector<unsigned int>, std::vector<unsigned int>>>& result)
{
    for (const auto& pair : result)
    {
        std::cout << "Left vector: ";
        for (const auto& leftValue : pair.first)
        {
            std::cout << leftValue << " ";
        }
        std::cout << std::endl;

        std::cout << "Right vector: ";
        for (const auto& rightValue : pair.second)
        {
            std::cout << rightValue << " ";
        }
        std::cout << std::endl;

        std::cout << std::endl;
    }
}
long getMemoryUse(){
    int who = RUSAGE_SELF;
    struct rusage usage;
    getrusage(who, &usage);
    return usage.ru_maxrss;
}
