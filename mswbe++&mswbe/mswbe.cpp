//
//  mswbe.cpp
//  
//
//  Created by 杨建业 on 2022/7/4.
//

#include "mswbe.h"
#include "Timer.h"


string Tools::integer_to_string(long long number) {
    std::vector<ui> sequence;
    if(number == 0) sequence.push_back(0);
    while(number > 0) {
        sequence.push_back(number%1000);
        number /= 1000;
    }
    
    char buf[5];
    std::string res;
    for(unsigned int i = sequence.size();i > 0;i --) {
        if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
        else sprintf(buf, ",%03u", sequence[i-1]);
        res += std::string(buf);
    }
    return res;
}

VertexDegree::VertexDegree(): vertex(-1), degree(-1){}

VertexDegree::VertexDegree(ui v, double d): vertex(v), degree(d){}

VertexDegree::~VertexDegree(){}

Node::Node(vector<ui> L, vector<ui> R, vector<ui> CR, vector<ui> XR, double minW, double maxW){
    this->L = L;
    this->R = R;
    this->CR = CR;
    this->XR = XR;
    this->minW = minW;
    this->maxW = maxW;
}

Node::~Node(){}

MSWBE::MSWBE(string path, double delta_value, int left_size, int right_size, int reduction_strategy, int sort_strategy){
    data_file_path = path;
    delta = delta_value;
    left_thd = left_size;
    right_thd = right_size;
    result_count = 0;
    graph_reduction_strategy = reduction_strategy;
    vertex_order_strategy = sort_strategy;
}

MSWBE::~MSWBE(){
    
}

void MSWBE::read_graph_txt(string file_path){
    cout<<"Start reading graph: "<<file_path<<endl;
    ifstream input_file(file_path, ios::in);
    map<ui, set<udp>> biG;
    if (!input_file.is_open()){
        cout << "Cannot open file : "<<file_path<<endl;exit(1);
    }
    else{
        input_file >> n1 >> n2 >> m;
        n = n1 + n2;
        cout<<"n1 = "<<n1<<", n2 = "<<n2<<", m = "<<m<<endl;
        ui u, v;
        double w;
        while (input_file >> u >> v >> w) {
            assert(u != v);
            assert(u >= 0 && u < n);
            assert(v >= 0 && v < n);
            assert(w > 0);
            biG[u].insert(make_pair(v, w));
            biG[v].insert(make_pair(u, w));
        }
        //cout << "n = " << n << endl;
        //cout << "biG.size() = " << biG.size() << endl;
        //assert(biG.size() == n);
        m = 0;
        for(auto e : biG) m += e.second.size();
        assert(m%2 == 0); m /= 2;
        input_file.close();
    }
    
    pstart = new ui[n+1];
    edges = new udp[2*m];
    degree = new int[n];
    TMPdeg = new int[n];
    twoDeg = new int[n];
    side = new int[n];
        
    pstart[0] = 0;
    for(ui i = 0; i < n; i++){
        const set<udp> & neighbors = biG[i];
        ui s_idx = pstart[i];
        for(auto e : neighbors) edges[s_idx++] = e;
        pstart[i+1] = s_idx;
        degree[i] = neighbors.size();
        TMPdeg[i] = neighbors.size();
        twoDeg[i] = 0;
    }
    assert(pstart[n] == 2*m);
    for(ui i = 0; i < n1; i++) side[i] = 0;
    for(ui i = n1; i < n; i++) side[i] = 1;
    
}


void MSWBE::core_based_reduction(){
    cout << "Start core-based reduction." << endl;
    //assert(left_thd >= 2 && right_thd >= 2);
    if(left_thd < 2 || right_thd < 2)
        return;
    
    int * removed = new int[n];
    memset(removed, 0, sizeof(int)*n);
    
    queue<ui> to_remove_vertices;
    for(ui i = 0; i < n; i++){          //cout << "TMPdeg[" << i << "]:" << TMPdeg[i] << endl;
        if((side[i] == 0 && TMPdeg[i] < right_thd) || (side[i] == 1 && TMPdeg[i] < left_thd)){
            to_remove_vertices.push(i);
            removed[i] = 1;
        }
    }
    int delcnt = 0;
    while(!to_remove_vertices.empty()){
        ui u = to_remove_vertices.front();      //cout << "delete : " << u << endl;
        to_remove_vertices.pop();
        del_ver[u] = 1;
        delcnt++;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i].first;
            if(del_ver[v] == 1 || removed[v] == 1)
                continue;
            --TMPdeg[v];
            if((side[v] == 0 && TMPdeg[v] == right_thd-1) || (side[v] == 1 && TMPdeg[v] == left_thd-1)){
                to_remove_vertices.push(v);
                removed[v] = 1;
            }
        }
    }
    
    cout<<"Delete "<< delcnt <<" vertices."<<endl;
    delete [] removed;
}


void MSWBE::weight_based_reduction(){
    cout << "Start weight-based reduction." << endl;
    //assert(left_thd >= 2 && right_thd >= 2);
    if(left_thd < 2 || right_thd < 2)
        return;
    
    vector< vector<udp> > tmpEdges(n);
    for(ui i = 0; i < n; i++){
        
        //cout << "u = " << i << ", ";
        
        for(ui j = pstart[i]; j < pstart[i+1]; j++){
            tmpEdges[i].emplace_back(edges[j]);
        }
        sort(tmpEdges[i].begin(), tmpEdges[i].end(), [](const udp &p1, const udp &p2) -> bool {
            return p1.second < p2.second || p1.second == p2.second && p1.first < p2.first;
        });
        
        //for(auto it = tmpEdges[i].begin(); it != tmpEdges[i].end(); it++){
        //    cout << it->first << ":" << it->second << ", ";
        //}
        //cout << endl;
    }
    
    int * removed = new int[n];
    memset(removed, 0, sizeof(int)*n);
    
    queue<ui> to_remove_vertices;
    for(ui i = 0; i < n; i++){          //cout << "TMPdeg[" << i << "]:" << TMPdeg[i] << endl;
        if(side[i] == 0 && !edgesSatisfied(tmpEdges[i], right_thd) || side[i] == 1 && !edgesSatisfied(tmpEdges[i], left_thd)){
            to_remove_vertices.push(i);
            removed[i] = 1;
        }
    }
    
    int delcnt = 0;
    while(!to_remove_vertices.empty()){
        ui u = to_remove_vertices.front();      //cout << "delete : " << u << endl;
        to_remove_vertices.pop();
        del_ver[u] = 1;
        delcnt++;
        for(ui i = pstart[u]; i < pstart[u+1]; i++){
            ui v = edges[i].first;
            if(del_ver[v] == 1 || removed[v] == 1)
                continue;
            --TMPdeg[v];
            
            // to remove edge <u,v> from edge list of v
            for(auto it = tmpEdges[v].begin(); it != tmpEdges[v].end(); it++){
                if(it->first == u){
                    tmpEdges[v].erase(it);
                    break;
                }
            }
            
            if(side[v] == 0 && !edgesSatisfied(tmpEdges[v], right_thd) || side[v] == 1 && !edgesSatisfied(tmpEdges[v], left_thd)){
                to_remove_vertices.push(v);
                removed[v] = 1;
            }
        }
    }
    
    cout<<"Delete "<< delcnt <<" vertices."<<endl;
    delete [] removed;
}

int MSWBE::edgesSatisfied(vector<udp>& P, int size_thd){
    int i=0, j=1;
    while(j < P.size()){
        if(P[j].second - P[i].second <= delta){
            j++;
        }
        else{
            if(j-i >= size_thd){
                return true;
            }
            do{
                i++;
            }
            while(P[j].second - P[i].second > delta);
        }
    }
    if(j-i >= size_thd){
        return true;
    }
    else{
        return false;
    }
}


vector<ui> MSWBE::get_two_hop_nei(ui v, vector<ui> &L, double lowW, double highW){
    //assert(v >=0 && v < r_n);
    vector<ui> two_hop_nei;
    ui *UV_cnt = new ui[r_n];
    memset(UV_cnt, 0, sizeof(ui)*r_n);
    vector<ui> may_2hop_list;
    
    for(auto u : L){
                                                        //cout << "u=" << u << endl;
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++){
            udp e = r_edges[i];                         //cout << "e : " << e.first << "," << e.second << endl;
            if(e.first == v)
                continue;
            if(e.second >= lowW && e.second <= highW) {
                if(UV_cnt[e.first] == 0){
                    may_2hop_list.emplace_back(e.first);
                }
                ++ UV_cnt[e.first];
            }
        }
    }
    
    for(auto e : may_2hop_list)
        if(UV_cnt[e] >= left_thd)
            two_hop_nei.emplace_back(e);
    delete [] UV_cnt;
    return two_hop_nei;
}


vector<ui> MSWBE::get_two_hop_nei(ui v){
    vector<ui> two_hop_nei;
    ui *UV_cnt = new ui[r_n];
    memset(UV_cnt, 0, sizeof(ui)*r_n);
    vector<ui> may_2hop_list;
    
    for(ui j = r_pstart[v]; j < r_pstart[v+1]; j++){
        ui u = r_edges[j].first;
        for(ui i = r_pstart[u]; i < r_pstart[u+1]; i++){
            udp e = r_edges[i];
            if(e.first == v)
                continue;
            if(UV_cnt[e.first] == 0) {
                may_2hop_list.emplace_back(e.first);
            }
            ++ UV_cnt[e.first];
        }
    }
    
    for(auto e : may_2hop_list){
        if(UV_cnt[e] >= left_thd) two_hop_nei.emplace_back(e);
    }
    delete [] UV_cnt;
    return two_hop_nei;
}

vector<ui> MSWBE::get_two_hop_nei(ui v, ui* &start_arr, udp* &edges_arr){
    vector<ui> two_hop_nei;
    ui *UV_cnt = new ui[r_n];
    memset(UV_cnt, 0, sizeof(ui)*r_n);
    vector<ui> may_2hop_list;
    
    for(ui j = start_arr[v]; j < start_arr[v+1]; j++){
        ui u = edges_arr[j].first;
        for(ui i = start_arr[u]; i < start_arr[u+1]; i++){
            udp e = edges_arr[i];
            if(e.first == v)
                continue;
            if(UV_cnt[e.first] == 0) {
                may_2hop_list.emplace_back(e.first);
            }
            ++ UV_cnt[e.first];
        }
    }
    
    for(auto e : may_2hop_list){
        if(UV_cnt[e] >= left_thd) two_hop_nei.emplace_back(e);
    }
    delete [] UV_cnt;
    return two_hop_nei;
}

vector<ui> MSWBE::get_two_hop_nei_with_range(ui v, ui* &start_arr, udp* &edges_arr, double &minW, double &maxW){
    
    vector<ui> two_hop_nei;
    ui *UV_cnt = new ui[r_n];
    memset(UV_cnt, 0, sizeof(ui)*r_n);
    vector<ui> may_2hop_list;
    
    double *min_arr = new double[r_n];
    double *max_arr = new double[r_n];
    for(ui i = 0; i < r_n; i++){
        min_arr[i] = INF;
        max_arr[i] = 0;
    }
    
    for(ui j = start_arr[v]; j < start_arr[v+1]; j++){
        ui u = edges_arr[j].first;
        for(ui i = start_arr[u]; i < start_arr[u+1]; i++){
            udp e = edges_arr[i];
            if(e.first == v)
                continue;
            if(UV_cnt[e.first] == 0) {
                may_2hop_list.emplace_back(e.first);
            }
            ++ UV_cnt[e.first];
            if(min_arr[e.first] > e.second){
                min_arr[e.first] = e.second;
            }
            if(max_arr[e.first] < e.second){
                max_arr[e.first] = e.second;
            }
        }
    }
    
    for(auto e : may_2hop_list){
        if(UV_cnt[e] >= left_thd){
            two_hop_nei.emplace_back(e);
            if(minW > min_arr[e]){
                minW = min_arr[e];
            }
            if(maxW < max_arr[e]){
                maxW = max_arr[e];
            }
        }
    }
    delete [] UV_cnt;
    delete [] min_arr;
    delete [] max_arr;
    return two_hop_nei;
}

void MSWBE::build_matrix_for_LCRXR(vector<ui> L, vector<ui> R, vector<ui> CR, vector<ui> XR){
#ifdef IS_DEBUGGING
    cout << "Enter matrix" << endl;
    cout << "L: "; print(L);
    cout << "R: "; print(R);
    cout << "CR: ";print(CR);
    cout << "XR: ";print(XR);
    cout << "Exit matrix" << endl;
#endif
    
    for(ui i = 0; i < max_CS_deg; i++)
        for(ui j = 0; j < max_CS_deg; j++)
            Matrix[i][j] = 0.0;
    int idx = 0;
    for(auto e : L) trans[e] = idx ++;
    for(auto e : R) trans[e] = idx ++;
    for(auto e : CR) trans[e] = idx ++;
    for(auto e : XR) trans[e] = idx ++;
    //L
    for(auto v : R) {
        for(ui i = r_pstart[v]; i < r_pstart[v+1]; i++) {
            ui u = r_edges[i].first;
            if(inL[u] == 1) {
                Matrix[trans[u]][trans[v]] = r_edges[i].second;
                Matrix[trans[v]][trans[u]] = r_edges[i].second;
            }
        }
    }
    for(auto v : CR) {
        for(ui i = r_pstart[v]; i < r_pstart[v+1]; i++) {
            ui u = r_edges[i].first;
            if(inL[u] == 1) {
                Matrix[trans[u]][trans[v]] = r_edges[i].second;
                Matrix[trans[v]][trans[u]] = r_edges[i].second;
            }
        }
    }
    for(auto v : XR) {
        for(ui i = r_pstart[v]; i < r_pstart[v+1]; i++) {
            ui u = r_edges[i].first;
            if(inL[u] == 1) {
                Matrix[trans[u]][trans[v]] = r_edges[i].second;
                Matrix[trans[v]][trans[u]] = r_edges[i].second;
            }
        }
    }
}


vector< pair<int, int> > MSWBE::find_sim_groups(vector< udp > &P, double d){
    assert(P.size() > 0);
    
    vector< pair<int, int> > ans;
    
    int i=0, j=1;
    while(j < P.size()){
        if(P[j].second - P[i].second <= delta){
            j++;
        }
        else{
            ans.emplace_back(make_pair(i, j-1));
            do{
                i++;
            }
            while(P[j].second - P[i].second > delta);
        }
    }
    ans.emplace_back(make_pair(i, j-1));
    return ans;
}


vector< vector<ui> > MSWBE::collect_branches(vector< pair<ui, ddp> > &P, double d){
    assert(P.size() > 0);
    
    vector< vector<ui> > ans;
    
    sort(P.begin(), P.end(), [](const pair<ui, ddp> &p1, const pair<ui, ddp> &p2) -> bool {
        return p1.second.first < p2.second.first || p1.second.first == p2.second.first && p1.second.second < p2.second.second;
    });
    
    for(int i = 0; i < P.size(); i++){
        vector<ui> L;
        L.emplace_back(P[i].first);
        for(int j = i+1; j < P.size(); j++){
            if(P[j].second.second - P[i].second.first <= d){
                L.emplace_back(P[j].first);
            }
            else if(P[j].second.first - P[i].second.first > d){
                break;
            }
        }
        
        sort(L.begin(), L.end());
        if(ans.empty()){
            ans.emplace_back(L);
        }
        else{
            vector<ui> L_prime = ans.back();
            int j,k;
            j = L_prime.size() - 1;
            k = L.size() - 1;
            int cnt = 0;
            while(j >= 0 && k >= 0){
                if(L_prime[j] == L[k]){
                    j--;
                    k--;
                    cnt++;
                }
                else if(L_prime[j] > L[k]){
                    j--;
                }
                else{
                    k--;
                    break;
                }
            }
            if(cnt < L.size()){
                ans.emplace_back(L);
            }
        }
        
    }
        
    return ans;
}

void MSWBE::sort_vertices(vector<ui> &ids, ui strategy){
    switch (strategy) {
        case 0: {// random vertex order
            cout << "Random vertex order is applied..." << endl;
            
            // we do nothing for random order
            
            break;
        }
        case 1: {// degree vertex order
            cout << "Degree vertex order is applied..." << endl;
            
            vector<VertexDegree> verdegs(ids.size());
            for(ui i = 0; i < ids.size(); i++){
                VertexDegree vd(ids[i], tmp_r_degree[ids[i]]);
                verdegs[i] = vd;
            }
            
            sort(verdegs.begin(), verdegs.end());
            
            for(ui i = 0; i < ids.size(); i++){
                ids[i] = verdegs[i].vertex;
            }
            
            break;
        }
        case 2:{// global summation based vertex order
            cout << "Global summation vertex order is applied..." << endl;
            vector<VertexDegree> verdegs(ids.size());
            for(ui i = 0; i < ids.size(); i++){
                ui v = ids[i];
                ui deg_sum = 0;
                for(ui j = tmp_r_pstart[v]; j < tmp_r_pstart[v+1]; j++){
                    ui u = tmp_r_edges[j].first;
                    deg_sum += tmp_r_degree[u];
                }
                VertexDegree vd(v, deg_sum);
                verdegs[i] = vd;
            }
            sort(verdegs.begin(), verdegs.end());
            
            for(ui i = 0; i < ids.size(); i++){
                ids[i] = verdegs[i].vertex;
            }
            
            break;
        }
        case 3:{// two hop neighbors based vertex order
            cout << "Two hop neighbors based order is applied..." << endl;
            vector<VertexDegree> verdegs(ids.size());
            for(ui i = 0; i < ids.size(); i++){
                ui v = ids[i];
                vector<ui> vec = get_two_hop_nei(v, tmp_r_pstart, tmp_r_edges);
                VertexDegree vd(v, vec.size());
                verdegs[i] = vd;
            }
            sort(verdegs.begin(), verdegs.end());
            
            for(ui i = 0; i < ids.size(); i++){
                ids[i] = verdegs[i].vertex;
            }
            break;
        }
        case 4:{// normalized two hop neighbors based vertex order
            cout << "Two hop neighbors based order is applied..." << endl;
            vector<VertexDegree> verdegs(ids.size());
            for(ui i = 0; i < ids.size(); i++){
                ui v = ids[i];
                double minW, maxW;
                minW = INF;
                maxW = 0;
                vector<ui> vec = get_two_hop_nei_with_range(v, tmp_r_pstart, tmp_r_edges, minW, maxW);
                double factor = (maxW - minW) / delta;
                double norm_d = (double)vec.size()/max(factor,1);
                VertexDegree vd(v, norm_d);
                verdegs[i] = vd;
            }
            sort(verdegs.begin(), verdegs.end());
            
            for(ui i = 0; i < ids.size(); i++){
                ids[i] = verdegs[i].vertex;
            }
            break;
        }
        default:{
            cout << "Please select the following vertex ordering: 0-random, 1-degree, 2-globalsum" << endl;
            return;
        }
    }
}


void MSWBE::advanced_mswbe(){
    cout << "advanced_mswbe is called." << endl;
    
    //vertex reduction
    del_ver = new bool[n];
    memset(del_ver, 0, sizeof(bool)*n);
    
    if(graph_reduction_strategy == 0){
        core_based_reduction();
    }
    else{
        weight_based_reduction();
    }
    
    //rebuild remaining graph
    r_n = 0; r_n1 = 0; r_n2 = 0;
    r_m = 0;
    for(ui u = 0; u < n1; u++) if(del_ver[u] == 0) {
        ++r_n; ++r_n1;
        r_m += TMPdeg[u];
    }
    for(ui u = n1; u < n; u++) if(del_ver[u] == 0) {
        ++r_n; ++r_n2;
        r_m += TMPdeg[u];
    }
    assert(r_m%2 == 0); r_m /= 2;
    
    cout<<"Reduced graph: n1 = "<<r_n1<<", n2 = "<<r_n2<<", m = "<<r_m<<endl;
    
    if(r_n == 0) return;
    
    
    tmp_oid = new ui[r_n];
    tmp_nid = new ui[n];
    
    ui vididx = 0;
    for(ui u = 0; u < n; u++)
        if(del_ver[u] == 0) {
            tmp_nid[u] = vididx;
            tmp_oid[vididx] = u;
            ++ vididx;
        }
    assert(r_n == vididx);
    
    tmp_r_pstart = new ui[r_n+1];
    tmp_r_edges = new udp[2*r_m];
    tmp_r_degree = new int[r_n];
    
    tmp_r_pstart[0] = 0;
    ui r_idx = 0;
    for(ui u = 0; u < r_n; u++){
        ui start_pos = tmp_r_pstart[r_idx];
        ui tdeg = 0;
        ui orivid = tmp_oid[u];
        for(ui i = pstart[orivid]; i < pstart[orivid+1]; i++){
            if(del_ver[edges[i].first] == 0) {
                tmp_r_edges[start_pos++] = make_pair(tmp_nid[edges[i].first], edges[i].second);
                ++ tdeg;
            }
        }
        tmp_r_degree[u] = tdeg;
        tmp_r_pstart[++r_idx] = start_pos;
    }
    
    assert(r_idx == r_n);

    // sort vertex
    vector<ui> vertices;
    for(ui v = (r_n1<=r_n2?0:r_n1); v < (r_n1<=r_n2?r_n1:r_n); v++){
        vertices.emplace_back(v);
    }

    sort_vertices(vertices, vertex_order_strategy);
    
    oid = new ui[r_n];
    nid = new ui[r_n];
    
    ui newid = r_n1<=r_n2?0:r_n1;
    for(ui v : vertices){
        nid[v] = newid;
        oid[newid] = v;
        ++ newid;
    }
    for(int i = r_n1<=r_n2?r_n1:0; i < (r_n1<=r_n2?r_n:r_n1); i++){
        oid[i] = i;
        nid[i] = i;
    }
    
    r_pstart = new ui[r_n+1];
    r_edges = new udp[2*r_m];
    r_degree = new int[r_n];
    
    r_pstart[0] = 0;
    r_idx = 0;
    for(ui u = 0; u < r_n; u++){
        ui start_pos = r_pstart[r_idx];
        ui tdeg = 0;
        ui orivid = oid[u];
        for(ui i = tmp_r_pstart[orivid]; i < tmp_r_pstart[orivid+1]; i++){
            r_edges[start_pos++] = make_pair(nid[tmp_r_edges[i].first], tmp_r_edges[i].second);
            ++ tdeg;
        }
        r_degree[u] = tdeg;
        r_pstart[++r_idx] = start_pos;
    }
    if(r_n1 <= r_n2){
        int tmpV = left_thd;
        left_thd = right_thd;
        right_thd = tmpV;
    }

    max_CS_deg = 0;

    for(ui v = (r_n1<=r_n2?0:r_n1); v < (r_n1<=r_n2?r_n1:r_n); v++){
        twoDeg[v] = get_two_hop_nei(v).size();
        if(twoDeg[v] + r_degree[v] > max_CS_deg)
            max_CS_deg = twoDeg[v] + r_degree[v];
    }
    
    max_CS_deg += 1;        //cout << "max_CS_deg = " << max_CS_deg << endl;
    
    Matrix = new double * [max_CS_deg];
    for(ui i = 0; i < max_CS_deg; i++)
        Matrix[i] = new double [max_CS_deg];
    trans = new ui[r_n];
    
    inL = new int[r_n];
    memset(inL, 0, sizeof(int)*r_n);
    deg_inL = new int[r_n];
    
    minmaxVec = new ddp[r_n];
    
    
    for(ui v = (r_n1<=r_n2?0:r_n1); v < (r_n1<=r_n2?r_n1:r_n); v++){
        
//#ifdef IS_DEBUGGING
        //cout << "dealing with v = " << v << endl;
//#endif
        
        vector<udp> P; // P stores all neighbors of v ignoring the weight constraint
        vector<ui> L, R, CR, XR;
        
        R.emplace_back(v);
        
        for(ui i = r_pstart[v]; i < r_pstart[v+1]; i++) {
            udp u = r_edges[i];
            P.emplace_back(u);
            L.emplace_back(u.first);
            inL[u.first] = 1;
            minmaxVec[u.first] = make_pair(u.second, u.second);
        }

        sort(P.begin(), P.end(), [](const udp &p1, const udp &p2) -> bool {
            return p1.second < p2.second || p1.second == p2.second && p1.first < p2.first;
        });
        
        
        if(L.size() < left_thd){
            for(auto e : L) inL[e] = 0;
            continue;
        }
        
        vector<ui> vec = get_two_hop_nei(v);    //print(vec);
        for(auto & w : vec){
            if(w > v){
                CR.emplace_back(w);
            }
            else if(w < v){
                XR.emplace_back(w);
            }
        }
        
        build_matrix_for_LCRXR(L, R, CR, XR);
        
        level_based_search(P, v);
        for(auto e : L) inL[e] = 0;
    }
}

void MSWBE::level_based_search(vector<udp> &P, ui vid){
    
    double minW, maxW;
#ifdef IS_DEBUGGING
    int level = 0;
#endif
    
    queue<Node> searchNodeQueue;
    
    vector< pair<int, int> > groupsInL = find_sim_groups(P, delta);
    
//#ifdef IS_DEBUGGING
    //cout << groupsInL.size() << endl;
//#endif
    
    vector<ui> R;
    R.emplace_back(vid);
    
    int max_level = 0;
    
    for(auto & g : groupsInL){
        vector<ui> L, CR, XR;
        minW = INF;
        maxW = 0;
        
                                //cout << "first = " << g.first << ", second=" << g.second << endl;
        
        for(int i = g.first; i <= g.second; i++){
            L.emplace_back(P[i].first);
        }
        minW = min(minW, P[g.first].second);
        maxW = max(maxW, P[g.second].second);
        
        if(L.size() < left_thd){
            continue;
        }
                                                                //print(L);
        
        double lowW = minW - delta;
        double highW = maxW + delta;
        
#ifdef IS_DEBUGGING
        cout << "min = " << minW << ", max = " << maxW << endl;
        cout << "low = " << lowW << ", high = " << highW << endl;
#endif
        
        vector<ui> vec = get_two_hop_nei(vid, L, lowW, highW);    //print(vec);
        sort(vec.begin(), vec.end());
        for(auto & w : vec){
            if(w > vid){
                CR.emplace_back(w);
            }
            else if(w < vid){
                XR.emplace_back(w);
            }
        }
        
        if(R.size() + CR.size() >= right_thd){
            Node node(L, R, CR, XR, minW, maxW);
            searchNodeQueue.emplace(node);
            max_level = max(max_level, R.size() + CR.size());
            
        }
    }
    
    vector< unordered_map<vector<ui>, vector<Node>, HashFunc, EqualFunc> > candidate_bicliques(max_level+1);
    
    //while(!searchNodeQueue.empty()){
    for(int level = 1; level <= max_level; level++){
        while(!searchNodeQueue.empty()){
            Node head = searchNodeQueue.front();
            searchNodeQueue.pop();
            
#ifdef IS_DEBUGGING
            cout << "---------------------new search node---------------------" << endl;
            print(head);
            cout << "---------------------------------------------------------" << endl;
#endif
            
            //if(head.L.size() >=  left_thd && head.R.size() >= right_thd){
                
                bool maximality = true;
                // if not exist v in CR or XR, such that LL \subseteq N(v),
                // we print the maximal biclique
                
                for(auto v : head.CR){           //cout << "CR : " << v << endl;
                    int cnt = 0;
                    double tmpMin = INF;
                    double tmpMax = 0;
                    for(auto u : head.L){
                        if(Matrix[trans[u]][trans[v]] > 0){
                            ++cnt;
                            tmpMin = min(tmpMin, Matrix[trans[u]][trans[v]]);
                            tmpMax = max(tmpMax, Matrix[trans[u]][trans[v]]);
                        }
                        else
                            break;
                    }
                    //if(cnt == LL.size() && (tmpMax - tmpMin <= delta)){
                    if(cnt == head.L.size() && (tmpMax - head.minW <= delta) && (head.maxW - tmpMin <= delta) && (tmpMax - tmpMin <= delta)){
                        maximality = false;
                        break;
                    }
                }
                for(auto v : head.XR){           //cout << "XR : " << v << endl;
                    int cnt = 0;
                    double tmpMin = INF;
                    double tmpMax = 0;
                    for(auto u : head.L){
                        if(Matrix[trans[u]][trans[v]] > 0){
                            ++cnt;
                            tmpMin = min(tmpMin, Matrix[trans[u]][trans[v]]);
                            tmpMax = max(tmpMax, Matrix[trans[u]][trans[v]]);
                        }
                        else
                            break;
                    }
                    
                    if(cnt == head.L.size() && (tmpMax - head.minW <= delta) && (head.maxW - tmpMin <= delta) && (tmpMax - tmpMin <= delta)){
                        maximality = false;
                        break;
                    }
                }
                if(maximality && head.R.size() >= right_thd)
                {
                    result_count++;
                    vector<ui>resL,resR;
                for(auto l:head.L){
                    resL.emplace_back(tmp_oid[oid[l]]);
                }//tmp_oid
                for(auto r:head.R){
                    resR.emplace_back(tmp_oid[oid[r]]);
                }
                sort(resL.begin(), resL.end());
                sort(resR.begin(), resR.end());

                if(r_n1 <= r_n2)
                    result_bicliques.emplace_back(make_pair(resR, resL));
                else
                    result_bicliques.emplace_back(make_pair(resL, resR));
#ifdef IS_DEBUGGING
                    cout << "-------------------collect 1 biclique---------------------" << endl;
                    print_biclique(head.L, head.R);
                    cout << "----------------------------------------------------------" << endl;
#endif
                    
#ifndef COUNT_ONLY
                    sort(head.L.begin(), head.L.end());
                    sort(head.R.begin(), head.R.end());
                    
                    if(r_n1 <= r_n2)
                        result_bicliques.emplace_back(make_pair(head.R, head.L));
                    else
                        result_bicliques.emplace_back(make_pair(head.L, head.R));
#endif
                }
            //}
            
            // step 1: for each v \in CR, expand R with v and branch on L \cap N(v), store candidate nodes into CB
            for(ui i = 0; i < head.CR.size(); i++){
                ui v = head.CR[i];
                
#ifdef IS_DEBUGGING
                cout << "cand v = " << v << endl;
#endif
                
                vector<ui> R_tmp;
                R_tmp = head.R;
                R_tmp.emplace_back(v);
                
                double lowW = head.minW - delta;
                double highW = head.maxW + delta;
                
                //cout << "lowW = " << lowW << ", highW = " << highW << endl;
                
                // get the neighbors of v in L with weight in range [low, high]
                vector< pair<ui, ddp> > P;
                
                for(auto x : head.L){
                    minmaxVec[x].first = INF;
                    minmaxVec[x].second = 0;
                    
                    ui neigCnt = 0;
                    for(auto y : R_tmp){
                        //cout << "x = " << x << ", y=" << y << ", Matrix = " << Matrix[trans[x]][trans[y]] << endl;
                        if(Matrix[trans[x]][trans[y]] > 0 && Matrix[trans[x]][trans[y]] >= lowW && Matrix[trans[x]][trans[y]] <= highW){
                            if(minmaxVec[x].first > Matrix[trans[x]][trans[y]])
                                minmaxVec[x].first = Matrix[trans[x]][trans[y]];
                            if(minmaxVec[x].second < Matrix[trans[x]][trans[y]])
                                minmaxVec[x].second = Matrix[trans[x]][trans[y]];
                            neigCnt++;
                        }
                    }
                    if(neigCnt >= R_tmp.size() && minmaxVec[x].second - minmaxVec[x].first <= delta)
                        P.emplace_back(make_pair(x, minmaxVec[x]));
                    
                }

                if(P.size() < left_thd){
                    continue;
                }
                
                vector< vector<ui> > commonNeighborSets = collect_branches(P, delta);
                
#ifdef IS_DEBUGGING
                cout << "Size of commonNeighborSets : " << commonNeighborSets.size() << endl;
                for(auto & L_prime : commonNeighborSets){
                    print(L_prime);
                }
#endif
                
                // branching the search by each of the groups in commonNeighborSets
                for(auto & L_prime : commonNeighborSets){
                    if(L_prime.size() < left_thd){
                        continue;
                    }
                    
#ifdef IS_DEBUGGING
                    cout << "branch on L: ";
                    print(L_prime);
#endif
                    vector<ui> R_prime, CR_prime, XR_prime;
                    R_prime = R_tmp;
                    
                    //vector<ui> XR_prime;
                    
                    minW = INF;
                    maxW = 0;
                    for(auto x : L_prime){
                        minW = min(minW, minmaxVec[x].first);
                        maxW = max(maxW, minmaxVec[x].second);
                    }
                    
                    lowW = minW - delta;
                    highW = maxW + delta;
                    
                    for(ui j = i+1; j < head.CR.size(); j++){
                        ui y = head.CR[j];
                        int da = 0;
                        int dc = 0;
                        for(auto z : L_prime){
                            if(Matrix[trans[z]][trans[y]] > 0 && Matrix[trans[z]][trans[y]] >= lowW && Matrix[trans[z]][trans[y]] <= highW){
                                ++dc;
                                //if(Matrix[trans[z]][trans[y]] >= minW && Matrix[trans[z]][trans[y]] <= maxW){
                                if(abs(Matrix[trans[z]][trans[v]]-Matrix[trans[z]][trans[y]]) <= EPSILON){
                                    ++da;
                                }
                            }
                        }
                        
                        if(dc >= left_thd){
                            if(da >= L_prime.size()){
                                R_prime.emplace_back(y);
                            }
                            else{
                                CR_prime.emplace_back(y);
                            }
                        }
                    }
                    
                    for(auto y : head.XR){
                        int d = 0;
                        for(auto z : L_prime){
                            if(Matrix[trans[z]][trans[y]] > 0 && Matrix[trans[z]][trans[y]] >= lowW && Matrix[trans[z]][trans[y]] <= highW){
                                ++d;
                            }
                        }
                        if(d >= left_thd){
                            XR_prime.emplace_back(y);
                        }
                    }
                    
                    if(R_prime.size() + CR_prime.size() >= right_thd){
                        Node newNode(L_prime, R_prime, CR_prime, XR_prime, minW, maxW);
                        
                        sort(R_prime.begin(), R_prime.end());
                        if(candidate_bicliques[R_prime.size()-1].find(R_prime) == candidate_bicliques[R_prime.size()-1].end()){
                            vector<Node> newList;
                            newList.emplace_back(newNode);
                            candidate_bicliques[R_prime.size()-1][R_prime] = newList;
                        }
                        else{
                            candidate_bicliques[R_prime.size()-1][R_prime].emplace_back(newNode);
                        }
                    }
                }
                
                // add v into XR
                head.XR.emplace_back(v);
                
#ifdef IS_DEBUGGING
                cout << "put " << v << " into xr" << endl;
#endif
            }
        }
        
#ifdef IS_DEBUGGING
        cout << "Finishing computing level : " << level << endl;
#endif
        
        // step 2: remove the subset in L side of bicliques in candidateBicliques and push into q1
        for(auto it = candidate_bicliques[level].begin(); it != candidate_bicliques[level].end(); it++){
            
            vector<Node> nodeList = it->second;
            sort(nodeList.begin(), nodeList.end(), [](const Node &p1, const Node &p2) -> bool {
                return p1.L.size() < p2.L.size();
            });
            
            ui i,j;
            for(i = 0; i < nodeList.size(); i++){
                
                vector<ui> Li = nodeList[i].L;
                for(j = i+1; j < nodeList.size(); j++){
                    vector<ui> Lj = nodeList[j].L;
                    int r,s;
                    r = s = 0;
                    int cnt = 0;
                    while(r < Li.size() && s < Lj.size()){
                        if(Li[r] == Lj[s]){
                            r++;
                            s++;
                            cnt++;
                        }
                        else if(Li[r] > Lj[s]){
                            s++;
                        }
                        else{
                            break;
                        }
                    }
                    if(cnt == Li.size()){
                        break;
                    }
                }
                if(j >= nodeList.size()){
                    searchNodeQueue.emplace(nodeList[i]);
                }
            }
        }
    }
}


void MSWBE::baseline_mswbe(){
    cout << "baseline_mswbe is called." << endl;
    
    //vertex reduction
    del_ver = new bool[n];
    memset(del_ver, 0, sizeof(bool)*n);
    
    if(graph_reduction_strategy == 0){
        core_based_reduction();
    }
    else{
        weight_based_reduction();
    }
    
    //rebuild remaining graph
    r_n = 0; r_n1 = 0; r_n2 = 0;
    r_m = 0;
    for(ui u = 0; u < n1; u++) if(del_ver[u] == 0) {
        ++r_n; ++r_n1;
        r_m += TMPdeg[u];
    }
    for(ui u = n1; u < n; u++) if(del_ver[u] == 0) {
        ++r_n; ++r_n2;
        r_m += TMPdeg[u];
    }
    assert(r_m%2 == 0); r_m /= 2;
    
    cout<<"Reduced graph: n1 = "<<r_n1<<", n2 = "<<r_n2<<", m = "<<r_m<<endl;
    
    if(r_n == 0) return;
    
    tmp_oid = new ui[r_n];
    tmp_nid = new ui[n];
    
    ui vididx = 0;
    for(ui u = 0; u < n; u++)
        if(del_ver[u] == 0) {
            tmp_nid[u] = vididx;
            tmp_oid[vididx] = u;
            ++ vididx;
        }
    assert(r_n == vididx);
    
    tmp_r_pstart = new ui[r_n+1];
    tmp_r_edges = new udp[2*r_m];
    tmp_r_degree = new int[r_n];
    
    tmp_r_pstart[0] = 0;
    ui r_idx = 0;
    for(ui u = 0; u < r_n; u++){
        ui start_pos = tmp_r_pstart[r_idx];
        ui tdeg = 0;
        ui orivid = tmp_oid[u];
        for(ui i = pstart[orivid]; i < pstart[orivid+1]; i++){
            if(del_ver[edges[i].first] == 0) {
                tmp_r_edges[start_pos++] = make_pair(tmp_nid[edges[i].first], edges[i].second);
                ++ tdeg;
            }
        }
        tmp_r_degree[u] = tdeg;
        tmp_r_pstart[++r_idx] = start_pos;
    }
    
    assert(r_idx == r_n);

    // sort vertex
    vector<ui> vertices;
    for(ui v = (r_n1<=r_n2?0:r_n1); v < (r_n1<=r_n2?r_n1:r_n); v++){
        vertices.emplace_back(v);
    }

    sort_vertices(vertices, vertex_order_strategy);
    
    oid = new ui[r_n];
    nid = new ui[r_n];
    
    ui newid = r_n1<=r_n2?0:r_n1;
    for(ui v : vertices){
        nid[v] = newid;
        oid[newid] = v;
        ++ newid;
    }
    for(int i = r_n1<=r_n2?r_n1:0; i < (r_n1<=r_n2?r_n:r_n1); i++){
        oid[i] = i;
        nid[i] = i;
    }
    
    r_pstart = new ui[r_n+1];
    r_edges = new udp[2*r_m];
    r_degree = new int[r_n];
    
    r_pstart[0] = 0;
    r_idx = 0;
    for(ui u = 0; u < r_n; u++){
        ui start_pos = r_pstart[r_idx];
        ui tdeg = 0;
        ui orivid = oid[u];
        for(ui i = tmp_r_pstart[orivid]; i < tmp_r_pstart[orivid+1]; i++){
            r_edges[start_pos++] = make_pair(nid[tmp_r_edges[i].first], tmp_r_edges[i].second);
            ++ tdeg;
        }
        r_degree[u] = tdeg;
        r_pstart[++r_idx] = start_pos;
    }
    
    if(r_n1 <= r_n2){
        int tmpV = left_thd;
        left_thd = right_thd;
        right_thd = tmpV;
    }
    
    max_CS_deg = 0;
    for(ui v = (r_n1<=r_n2?0:r_n1); v < (r_n1<=r_n2?r_n1:r_n); v++){

        twoDeg[v] = get_two_hop_nei(v).size();
        if(twoDeg[v] + r_degree[v] > max_CS_deg)
            max_CS_deg = twoDeg[v] + r_degree[v];
    }
    
    max_CS_deg += 1;        //cout << "max_CS_deg = " << max_CS_deg << endl;
    
    Matrix = new double * [max_CS_deg ];
    for(ui i = 0; i < max_CS_deg; i++)
        Matrix[i] = new double [max_CS_deg ];
    trans = new ui[r_n];
    
    inL = new int[r_n];
    memset(inL, 0, sizeof(int)*r_n);
    deg_inL = new int[r_n];
    
    minmaxVec = new ddp[r_n];
    
    //cout<<"\tStart enumeration procedure: "<<endl;
    // we enumerate by anchoring on the right side
    
    for(ui v = (r_n1<=r_n2?0:r_n1); v < (r_n1<=r_n2?r_n1:r_n); v++){
        
//#ifdef IS_DEBUGGING
        //cout << "dealing with v = " << v << endl;
//#endif
        
        //find maximal similar-weight bicliques containing u
        vector<ui> L, R, CR, XR;
        R.emplace_back(v);
        
        for(ui i = r_pstart[v]; i < r_pstart[v+1]; i++){
            udp u = r_edges[i];
            L.emplace_back(u.first);
            inL[u.first] = 1;
        }
        
        if(L.size() < left_thd){
            for(auto e : L) inL[e] = 0;
            continue;
        }
        
        vector<ui> vec = get_two_hop_nei(v);    //print(vec);
        for(auto & w : vec){
            if(w > v){
                CR.emplace_back(w);
            }
            else if(w < v){
                XR.emplace_back(w);
            }
        }
        
        build_matrix_for_LCRXR(L, R, CR, XR);
        
        branch_and_bound(L, R, CR, XR);
        for(auto e : L) inL[e] = 0;
        
    }
}

void MSWBE::branch_and_bound(vector<ui> &L, vector<ui> &R, vector<ui> &CR, vector<ui> &XR){
#ifdef IS_DEBUGGING
    cout << "L: "; print(L);
    cout << "R: "; print(R);
    cout << "CR: ";print(CR);
    cout << "XR: ";print(XR);
#endif
    
    
    // Step 1: spliting L by weight constraints, and collecting results
    if(L.size() >=  left_thd && R.size() >= right_thd){
        vector< pair<ui, ddp> > P;
        
        for(auto x : L){
            minmaxVec[x].first = INF;
            minmaxVec[x].second = 0;
            //ui neigCnt = 0;
            for(auto y : R){
                //cout << "x = " << x << ", y=" << y << ", Matrix = " << Matrix[trans[x]][trans[y]] << endl;
                if(Matrix[trans[x]][trans[y]] > 0){
                    if(minmaxVec[x].first > Matrix[trans[x]][trans[y]])
                        minmaxVec[x].first = Matrix[trans[x]][trans[y]];
                    if(minmaxVec[x].second < Matrix[trans[x]][trans[y]])
                        minmaxVec[x].second = Matrix[trans[x]][trans[y]];
                }
            }
            //cout << "x = " << x << ": " << minmaxVec[x].first << "," << minmaxVec[x].second << endl;
            //if(neigCnt >= right_thd && minmaxVec[x].second - minmaxVec[x].first <= delta)
            if(minmaxVec[x].second - minmaxVec[x].first <= delta)
                P.emplace_back(make_pair(x, minmaxVec[x]));
        }
        
        if(P.size() >= left_thd){
            
            vector< vector<ui> > commonNeighborSets = collect_branches(P, delta);
            
            #ifdef IS_DEBUGGING
            cout << "Size of commonNeighborSets : " << commonNeighborSets.size() << endl;
            for(auto & LL : commonNeighborSets){
                print(LL);
            }
            #endif
            
            for(auto & LL : commonNeighborSets){
                if(LL.size() < left_thd){
                    continue;
                }
                
                double minW = INF;
                double maxW = 0;
                for(auto x : LL){
                    minW = min(minW, minmaxVec[x].first);
                    maxW = max(maxW, minmaxVec[x].second);
                }
                bool maximality = true;
                // if not exist v in CR or XR, such that LL \subseteq N(v),
                // we print the maximal biclique
                
                for(auto v : CR){           //cout << "CR : " << v << endl;
                    int cnt = 0;
                    double tmpMin = INF;
                    double tmpMax = 0;
                    for(auto u : LL){
                        if(Matrix[trans[u]][trans[v]] > 0){
                            ++cnt;
                            tmpMin = min(tmpMin, Matrix[trans[u]][trans[v]]);
                            tmpMax = max(tmpMax, Matrix[trans[u]][trans[v]]);
                        }
                        else
                            break;
                    }
                    //if(cnt == LL.size() && (tmpMax - tmpMin <= delta)){
                    if(cnt == LL.size() && (tmpMax - minW <= delta) && (maxW - tmpMin <= delta)&& (tmpMax - tmpMin <= delta)){
                        maximality = false;
                        break;
                    }
                }
                for(auto v : XR){           //cout << "XR : " << v << endl;
                    int cnt = 0;
                    double tmpMin = INF;
                    double tmpMax = 0;
                    for(auto u : LL){
                        if(Matrix[trans[u]][trans[v]] > 0){
                            ++cnt;
                            tmpMin = min(tmpMin, Matrix[trans[u]][trans[v]]);
                            tmpMax = max(tmpMax, Matrix[trans[u]][trans[v]]);
                        }
                        else
                            break;
                    }
                    //if(cnt == LL.size() && (tmpMax - tmpMin <= delta)){
                    if(cnt == LL.size() && (tmpMax - minW <= delta) && (maxW - tmpMin <= delta) && (tmpMax - tmpMin <= delta)){
                        maximality = false;
                        break;
                    }
                }
                if(maximality && R.size() >= right_thd)
                {
                    result_count++;
                    
                    #ifdef IS_DEBUGGING
                    cout << "-------------------collect 1 biclique---------------------" << endl;
                    print_biclique(LL,R);
                    cout << "----------------------------------------------------------" << endl;
                    #endif
                    
                    #ifndef COUNT_ONLY
                    sort(LL.begin(), LL.end());
                    sort(R.begin(), R.end());
                        
                    if(r_n1 <= r_n2)
                        result_bicliques.emplace_back(make_pair(R, LL));
                    else
                        result_bicliques.emplace_back(make_pair(LL, R));
                        
                    #endif
                }
            }
        }
        
    }
    
    // Step 2: Branching on CR
    for(ui i = 0; i < CR.size(); i++){
        ui v = CR[i];
        
        #ifdef IS_DEBUGGING
        cout << "cand v = " << v << endl;
        #endif
        
        vector<ui> L_prime, R_prime, CR_prime, XR_prime;
        
        R_prime = R;
        R_prime.emplace_back(v);
        
        for(auto u : L){
            if(Matrix[trans[u]][trans[v]] > 0){
                L_prime.emplace_back(u);
            }
        }
        
        
        /*for(ui j = i+1; j < CR.size(); j++){
            ui y = CR[j];
            ui da = 0;
            ui dc = 0;
            for(auto z : L_prime){
                if(Matrix[trans[z]][trans[y]] > 0){
                    ++dc;
                    //if(abs(Matrix[trans[z]][trans[v]]-Matrix[trans[z]][trans[y]]) <= EPSILON){
                    //    ++da;
                    //}
                }
            }
            if(dc >= left_thd){
                //if(da >= L_prime.size()){
                //    R_prime.emplace_back(y);
                //}
                //else{
                    CR_prime.emplace_back(y);
                //}
            }
        }*/
        
        
        
        for(ui j = i+1; j < CR.size(); j++){
            ui y = CR[j];
            ui da = 0;
            ui dc = 0;
            for(auto z : L_prime){
                if(Matrix[trans[z]][trans[y]] > 0){
                    ++dc;
                    if(abs(Matrix[trans[z]][trans[v]] - Matrix[trans[z]][trans[y]]) <= EPSILON){
                        ++da;
                    }
                }
            }
            if(dc >= left_thd){
                if(da >= L_prime.size()){
                    R_prime.emplace_back(y);
                }
                else{
                    CR_prime.emplace_back(y);
                }
            }
        }
        
        for(auto y : XR){
            int d = 0;
            for(auto z : L_prime){
                if(Matrix[trans[z]][trans[y]] > 0){
                    ++d;
                }
            }
            if(d >= left_thd){
                XR_prime.emplace_back(y);
            }
        }
        
        if(R_prime.size() + CR_prime.size() >= right_thd && L_prime.size() >= left_thd){
            branch_and_bound(L_prime, R_prime, CR_prime, XR_prime);
        }
        
        XR.emplace_back(v);
        
#ifdef IS_DEBUGGING
        cout << "put " << v <<" into xr" << endl;
#endif
        
    }
    
}

long MSWBE::getMemoryUse(){
    int who = RUSAGE_SELF;
    struct rusage usage;
    getrusage(who, &usage);
    return usage.ru_maxrss;
}

void MSWBE::delete_memory(){
    if(pstart != nullptr){
        delete [] pstart;
        pstart = nullptr;
    }
    if(edges != nullptr){
        delete [] edges;
        edges = nullptr;
    }
    if(degree != nullptr){
        delete [] degree;
        degree = nullptr;
    }
    if(TMPdeg != nullptr){
        delete [] TMPdeg;
        TMPdeg = nullptr;
    }
    if(twoDeg != nullptr){
        delete [] twoDeg;
        twoDeg = nullptr;
    }
    if(side != nullptr){
        delete [] side;
        side = nullptr;
    }
    if(r_pstart != nullptr){
        delete [] r_pstart;
        r_pstart = nullptr;
    }
    if(r_edges != nullptr){
        delete [] r_edges;
        r_edges = nullptr;
    }
    if(r_degree != nullptr){
        delete [] r_degree;
        r_degree = nullptr;
    }
    if(oid != nullptr){
        delete [] oid;
        oid = nullptr;
    }
    if(nid != nullptr){
        delete [] nid;
        nid = nullptr;
    }
    if(tmp_r_pstart != nullptr){
        delete [] tmp_r_pstart;
        tmp_r_pstart = nullptr;
    }
    if(tmp_r_edges != nullptr){
        delete [] tmp_r_edges;
        tmp_r_edges = nullptr;
    }
    if(tmp_r_degree != nullptr){
        delete [] tmp_r_degree;
        tmp_r_degree = nullptr;
    }
    if(tmp_oid != nullptr){
        delete [] tmp_oid;
        tmp_oid = nullptr;
    }
    if(tmp_nid != nullptr){
        delete [] tmp_nid;
        tmp_nid = nullptr;
    }
    for(int i = 0; i < max_CS_deg; i++) if(Matrix[i] != nullptr){
        delete [] Matrix[i];
        Matrix[i] = nullptr;
    }
    if(trans != nullptr){
        delete [] trans;
        trans = nullptr;
    }
    if(del_ver != nullptr){
        delete [] del_ver;
        del_ver = nullptr;
    }
    if(inL != nullptr){
        delete [] inL;
        inL = nullptr;
    }
    if(deg_inL != nullptr){
        delete [] deg_inL;
        deg_inL = nullptr;
    }
}



void MSWBE::print(vector<ui> &vec){
    for(auto e : vec){
        cout << e << " ";
    }
    cout << endl;
}

void MSWBE::print_biclique(vector<ui> &L, vector<ui> &R){
    for(auto e : L){
        cout << e << " ";
    }
    cout << "| ";
    for(auto e : R){
        cout << e << " ";
    }
    cout << endl;
}
        
void MSWBE::print_biclique(ofstream &fout, vector<ui> &L, vector<ui> &R){
    for(auto e : L){
        fout << e << " ";
    }
    fout << "| ";
    for(auto e : R){
        fout << e << " ";
    }
    fout << endl;
}

void MSWBE::print(vector< pair<ui, ddp> > &P){
    for(auto e : P){
        cout << e.first << ": " << e.second.first << "-" << e.second.second << endl;
    }
    cout << endl;
}

void MSWBE::print_results(){
    cout << "\tTotal # results : " << result_count << endl;
    #ifndef COUNT_ONLY
    for(auto e : result_bicliques){
        print_biclique(e.first, e.second);
    }
    #endif
}

void MSWBE::check_results(){
    ofstream fout("../experiment/output.txt");
    fout << "\tTotal # results : " << result_count << endl;
    for(int i = 0; i < result_bicliques.size(); i++){
        sort(result_bicliques[i].first.begin(), result_bicliques[i].first.end());
        sort(result_bicliques[i].second.begin(), result_bicliques[i].second.end());
        #ifndef COUNT_ONLY
            print_biclique(fout, result_bicliques[i].first, result_bicliques[i].second);
        #endif
    }

    for(int i = 0; i < result_bicliques.size(); i++){
        for(int j = 0; j < result_bicliques.size(); j++){
            if(i != j && isSubsetOf(result_bicliques[i].first, result_bicliques[j].first) && isSubsetOf(result_bicliques[i].second, result_bicliques[j].second)){
                cout << i << ", " << j << endl;
                print_biclique(result_bicliques[i].first, result_bicliques[i].second);
                print_biclique(result_bicliques[j].first, result_bicliques[j].second);
            }
        }
    }
}

bool MSWBE::isSubsetOf(vector<ui> &Sub, vector<ui> &Sup){
    int r,s;
    r = s = 0;
    
    while(r < Sub.size() && s < Sup.size()){
        if(Sub[r] == Sup[s]){
            r++;
            s++;
        }
        else if(Sub[r] < Sup[s]){
            return false;
        }
        else{
            s++;
        }
    }
    if(r == Sub.size()){
        return true;
    }else{
        return false;
    }
}
        
bool MSWBE::isSameBiclique(pair<vector<ui>, vector<ui> > &b1, pair<vector<ui>, vector<ui> > &b2){
    return isSameVector(b1.first, b2.first) && isSameVector(b1.second, b2.second);
}
        
bool MSWBE::isSameVector(vector<ui> &vec1, vector<ui> &vec2){
    if(vec1.size() != vec2.size()){
        return false;
    }
    for(int i=0; i<vec1.size(); i++){
        if(vec1[i] != vec2[i]){
            return false;
        }
    }
    return true;
}

void MSWBE::print(vector<udp> &vec){
    for(auto e : vec){
        cout << e.first << ":" << e.second << " ";
    }
    cout << endl;
}

void MSWBE::print(Node &node){
    cout << "L: "; print(node.L);
    cout << "R: "; print(node.R);
    cout << "CR: ";print(node.CR);
    cout << "XR: ";print(node.XR);
}

void printResultBicliques(const std::vector<std::pair<std::vector<ui>, std::vector<ui>>>& result_bicliques) {
    // 遍历 result_bicliques 中的每个元素
    for (const auto& biclique : result_bicliques) {
        // 打印第一个向量
        std::cout << "First vector: ";
        for (const auto& val : biclique.first) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // 打印第二个向量
        std::cout << "Second vector: ";
        for (const auto& val : biclique.second) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << std::endl;
    }
}
void run_mswbe(int argc, char *argv[]){
    
    MSWBE* mswbe = NULL;
    cout << "argc = " << argc << endl;
    if(argc < 5){
        cout << "Too few arguments" << endl;
        return;
    }
    else if(argc == 5){
        for( int i = 0; i < argc; ++i )
            printf( "argv[%d]=%s\n", i, argv[i] );
        mswbe = new MSWBE(argv[1], atof(argv[2]), atoi(argv[3]), atoi(argv[4]));
    }
    else if(argc == 6){
        for( int i = 0; i < argc; ++i )
            printf( "argv[%d]=%s\n", i, argv[i] );
        mswbe = new MSWBE(argv[1], atof(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    }
    else{
        for( int i = 0; i < argc; ++i )
            printf( "argv[%d]=%s\n", i, argv[i] );
        mswbe = new MSWBE(argv[1], atof(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
    }
    
    mswbe->read_graph_txt(argv[1]);
    
    cout<<"delta="<<argv[2]<<", "<<"left_thd="<<argv[3]<<", "<<"right_thd="<<argv[4]<<endl;
    
    Timer t;
    if(argc >= 8 && atoi(argv[7]) == 0){
        mswbe->baseline_mswbe();
    }
    else{
        mswbe->advanced_mswbe();
    }
    //mswbe.result_bicliques.clear();
    
    cout<<"\t------------------------------------------------"<<endl;
    //cout<<"\tNumber of Similar-weight Bicliques: "<<mswbe.result_count<<endl;
    mswbe->print_results();
    cout<<"\tTime cost (without I/O): "<<Tools::integer_to_string(t.elapsed())<<"(s) or "<<Tools::integer_to_string(t.elapsed_in_millisec())<<"(ms)"<<endl;
    cout << "\tMemory usage : " << mswbe->getMemoryUse() / (1024) << "(MB)" << endl;
    cout<<"\t------------------------------------------------"<<endl;
    printResultBicliques(mswbe->result_bicliques);
    
    //mswbe.check_results();
    
    mswbe->delete_memory();
    //return mswbe->result_bicliques;
}
        
int main(int argc, char *argv[]){
    
    run_mswbe(argc, argv);
    
    //vector< pair<vector<ui>, vector<ui> > > baseline_result_bicliques = run_mswbe(argv);
    //vector< pair<vector<ui>, vector<ui> > > advanced_result_bicliques = run_mswbe(argv);
    
    /*for(int i = 0; i < baseline_result_bicliques.size(); i++){
        int j = 0;
        for(; j < advanced_result_bicliques.size(); j++){
            if(MSWBE::isSameBiclique(advanced_result_bicliques[j], baseline_result_bicliques[i]))
                break;
        }
        if(j >= advanced_result_bicliques.size()){
            MSWBE::print_biclique(baseline_result_bicliques[i].first, baseline_result_bicliques[i].second);
        }
    }
    cout << baseline_result_bicliques.size() << ", " << advanced_result_bicliques.size() << endl;
    */

    return 0;
}
