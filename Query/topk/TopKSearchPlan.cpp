
#include "TopKSearchPlan.h"

using namespace std;
RootStrategy global_root_strategy = RootStrategy::MinCand;
string root_var_name = "";
bool global_root_strategy_enable = false;
string RootStrategyToString(RootStrategy strategy)
{
  string r = "Invalid Strategy";
  switch (strategy) {
    case RootStrategy::Random:
      r = "Random";
      break;
    case RootStrategy::MinCand:
      r = "MinCand";
      break;
    case RootStrategy::MinDepth:
      r = "MinDepth";
      break;
    case RootStrategy::Ours:
      r = "Ours";
      break;
  }
  return r;

}

/**
 * count the tree depth regarding the root_id as root.
 * @note must be a tree, not a graph with loops
 * @param neighbours the tree structure
 * @param root_id root
 * @param total_vars_num the node number of the tree
 * @return the max tree depth
 */
std::size_t TopKSearchPlan::CountDepth(map<TYPE_ENTITY_LITERAL_ID, vector<TYPE_ENTITY_LITERAL_ID>>& neighbours,
                                       TYPE_ENTITY_LITERAL_ID root_id,
                                       std::size_t total_vars_num,
                                       set<int> &property_ids) {
  auto vars_used_vec = new bool[total_vars_num];
  auto vars_depth = new std::size_t [total_vars_num];

  for(int i = 0;i<total_vars_num;i++)
    vars_used_vec[i] = false;
  for(int i = 0;i<total_vars_num;i++)
    vars_depth[i] = 0;

  vars_used_vec[root_id] = true;
  vars_depth[root_id]= 1;
  size_t max_depth = 1;
  std::stack<int> explore_id;
  explore_id.push(root_id);
//  cerr<<"-------------"<<endl;
//  cerr<<"depth["<<root_id<<"]="<<1<<endl;
  while(!explore_id.empty())
  {
    auto now_id = explore_id.top();
    explore_id.pop();
    for(const auto& child_id: neighbours[now_id])
    {
      if(vars_used_vec[child_id])
        continue;
      vars_used_vec[child_id] = true;
      if(property_ids.find(child_id)!=property_ids.end())
        continue;

      vars_depth[child_id] = vars_depth[now_id] + 1;
//      cerr<<now_id<<"->"<<child_id<<" "<<vars_depth[now_id] + 1<<endl;
      max_depth = std::max(max_depth,vars_depth[child_id]);
//      cerr<<"depth["<<child_id<<"]="<<vars_depth[child_id]<<endl;
      explore_id.push(child_id);
    }
  }
//  cerr<<"-------------"<<endl;
  delete[] vars_used_vec;
  delete[] vars_depth;

  return max_depth;
}

/**
 * count the tree depth regarding the root_id as root.
 * @note must be a tree, not a graph with loops
 * @param neighbours the tree structure
 * @param root_id root
 * @param total_vars_num the node number of the tree
 * @return the min distance from root node to leaves
 */
std::size_t TopKSearchPlan::CountDepthMin(map<TYPE_ENTITY_LITERAL_ID, vector<TYPE_ENTITY_LITERAL_ID>>& neighbours,
                                       TYPE_ENTITY_LITERAL_ID root_id,
                                       std::size_t total_vars_num) {
  auto vars_used_vec = new bool[total_vars_num];
  auto vars_depth = new std::size_t [total_vars_num];
  memset(vars_used_vec,0,total_vars_num);
  memset(vars_depth,0,total_vars_num);
  vars_used_vec[root_id] = true;
  vars_depth[root_id]= 1;
  std::stack<int> explore_id;
  explore_id.push(root_id);
  size_t min_depth = 1 << 30;
  while(!explore_id.empty())
  {
    auto now_id = explore_id.top();
    bool leaf_node = true;
    for(const auto& child_id: neighbours[now_id])
    {
      if(vars_used_vec[child_id])
        continue;
      vars_used_vec[child_id] = true;
      leaf_node = false;
      vars_depth[child_id] = vars_depth[now_id] + 1;
      explore_id.push(child_id);
    }
    if(leaf_node)
      min_depth = std::min(min_depth,vars_depth[now_id]);
    explore_id.pop();
  }
  delete[] vars_used_vec;
  delete[] vars_depth;

  return min_depth;
}



TopKSearchPlan::TopKSearchPlan(shared_ptr<BGPQuery> bgp_query, KVstore *kv_store,
                               Statistics *statistics, const QueryTree::Order& expression,
                               IDCachesSharePtr id_caches):
                               tree_root_(nullptr)
{
  this->is_cycle_graph_ = false;
  this->total_vars_num_ = bgp_query->get_total_var_num();
  TYPE_ENTITY_LITERAL_ID CONSTANT = -1;
  // a neighbour may have more than one edge
  // constructing the information structure
  for(decltype(total_vars_num_) i = 0; i< total_vars_num_ ; i++)
  {
    map<TYPE_ENTITY_LITERAL_ID,size_t> neighbour_position;

    neighbours_[i] = decltype(neighbours_.begin()->second)();
    predicates_constant_[i] = decltype(predicates_constant_.begin()->second)();
    predicates_ids_[i] = decltype(predicates_ids_.begin()->second)();
    directions_[i] = decltype(directions_.begin()->second)();

    auto var_descriptor = bgp_query->get_vardescrip_by_index(i);
    if(var_descriptor->var_type_ == VarDescriptor::VarType::Predicate)
      continue;

    auto &edge_index_list = var_descriptor->so_edge_index_;

    auto v_name = bgp_query->get_var_name_by_id(i);
    for(auto edge_id:edge_index_list) {
      decltype(i) nei_id;
      auto &triple = bgp_query->get_triple_by_index(edge_id);
      decltype(triple.object) nei_name;
      auto direction= TopKPlanUtil::EdgeDirection::NoEdge;
      if(triple.subject==v_name) {
        nei_name = triple.object;
        direction = TopKPlanUtil::EdgeDirection::OUT;
      }
      else {
        nei_name = triple.subject;
        direction = TopKPlanUtil::EdgeDirection::IN;
      }
      if(nei_name[0]!='?')
        nei_id = CONSTANT;
      else
        nei_id = bgp_query->get_var_id_by_name(nei_name);

      if(nei_id==CONSTANT)
        continue;

      TYPE_ENTITY_LITERAL_ID predicates_id = CONSTANT;
      auto predicate_constant = triple.predicate[0]!='?';
      if(predicate_constant)
        predicates_id = kv_store->getIDByPredicate(triple.predicate);
      else
        predicates_id = bgp_query->get_var_id_by_name(triple.predicate);

      size_t neighbour_pos;
      if(neighbour_position.find(nei_id)==neighbour_position.end())
      {
        neighbours_[i].push_back(nei_id);
        neighbour_pos = neighbour_position.size();
        neighbour_position[nei_id] = neighbour_pos;

        predicates_constant_[i].emplace_back();
        predicates_ids_[i].emplace_back();
        directions_[i].emplace_back();
      }


      predicates_constant_[i][neighbour_pos].push_back(predicate_constant);
      predicates_ids_[i][neighbour_pos].push_back(predicates_id);
      directions_[i][neighbour_pos].push_back(direction);
    }
  }
}


void TopKSearchPlan::GetPlan(shared_ptr<BGPQuery> bgp_query,
                             KVstore *kv_store,
                             Statistics *statistics,
                             const QueryTree::Order& expression,
                             IDCachesSharePtr id_caches,
                             std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients) {
  auto first_cycle = this->FindCycle();
  if (!first_cycle.empty())
    this->CutCycle(bgp_query,kv_store,statistics,id_caches);

  set<int> property_ids;
  for (auto& var_co_pair : *var_coefficients)
    if (var_co_pair.second != 0)
      property_ids.insert(bgp_query->get_var_id_by_name(var_co_pair.first));

  int selected_root = -1;
  if (!global_root_strategy_enable)
    selected_root = SelectRoot(bgp_query, id_caches,property_ids);
  else {
    switch (global_root_strategy) {
      case RootStrategy::Ours:
        selected_root = SelectRoot(bgp_query, id_caches,property_ids);
//        selected_root = this->SelectRootCostModelOurs(bgp_query, id_caches, total_vars_num_, k,property_ids);
        break;
      case RootStrategy::MinCand:
        selected_root = this->SelectRootMinCand(bgp_query,id_caches,property_ids);
        break;
      case RootStrategy::MinDepth:
        selected_root = this->SelectRootMinDepth(bgp_query,id_caches,property_ids);
        break;
      case RootStrategy::Random:
        selected_root = this->SelectRootRandom(bgp_query,id_caches,property_ids);
        break;
    }
  }
//  cout << RootStrategyToString(global_root_strategy) << "selected_root: " << selected_root << "--------------" << endl;

#ifdef TOPK_DEBUG_INFO
  cout << "chosen root node:" << bgp_query->get_var_name_by_id(min_score_root)<<"["<<min_score_root<<"]" << endl;
#endif
  BuildPlanFromRoot(selected_root);
}

int TopKSearchPlan::SelectRoot(shared_ptr<BGPQuery> &bgp_query, IDCachesSharePtr &id_caches,set<int> &property_ids) {
  int min_score_root = -1;
  auto min_score = DBL_MAX;

  // get the root id by arg min( candidates * tree depth )
  // if all nodes don't have candidates, choose the min tree depth
  for(unsigned i = 0; i < total_vars_num_; ++i) {
    auto var_descriptor = bgp_query->get_vardescrip_by_index(i);
    if(var_descriptor->var_type_ == VarDescriptor::VarType::Predicate)
      continue;

    if(property_ids.find(i) != property_ids.end())
      continue;

    if(id_caches->size() && id_caches->find(i) == id_caches->end())
      continue;

    // count the depth if regard node i as root
    auto tree_depth = CountDepth(neighbours_, i, total_vars_num_,property_ids);
    auto candidates_size = id_caches->find(i)->second->size();

    double score;
    if(id_caches->size() > 0) score = tree_depth * candidates_size;
    else score = tree_depth;  // in fact, empty candidate set means that no results for this query

    if (score < min_score) {
      min_score_root = i;
      min_score = score;
    }
  }
  return min_score_root;
}

int TopKSearchPlan::SelectRootMinDepth(shared_ptr<BGPQuery> &bgp_query,
                                       IDCachesSharePtr &id_caches,
                                       set<int> &property_ids) {
  int min_score_root = -1;
  int min_depth = 1<<20;

  // get the root id by arg min( candidates * tree depth )
  // if all nodes don't have candidates, choose the min tree depth
  for(decltype(total_vars_num_)  i = 0; i< total_vars_num_; i++)
  {
    auto var_descriptor = bgp_query->get_vardescrip_by_index(i);
    if(var_descriptor->var_type_ == VarDescriptor::VarType::Predicate)
      continue;

    if(property_ids.find(i)!=property_ids.end()) continue;

    if(id_caches->size()&&id_caches->find(i)==id_caches->end())
      continue;
// count the its depth if regard node i as root
    auto tree_depth = CountDepth(neighbours_, i, total_vars_num_,property_ids);
    //cerr<<"Min Depth["<<i<<"] "<< bgp_query->get_var_name_by_id(i)<<" depth="<<tree_depth<<endl;
    if(min_score_root==-1 || tree_depth <= min_depth)
    {
      min_score_root = i;
      min_depth = tree_depth;
    }
  }
  return min_score_root;
}

int TopKSearchPlan::SelectRootMinCand(shared_ptr<BGPQuery> &bgp_query,
                                       IDCachesSharePtr &id_caches,
                                       set<int> &property_ids) {
  int min_cand_root = -1;
  int min_cand = INT_MAX;

  for(decltype(total_vars_num_)  i = 0; i< total_vars_num_; i++) {
    auto var_descriptor = bgp_query->get_vardescrip_by_index(i);
    if (var_descriptor->var_type_ == VarDescriptor::VarType::Predicate)
      continue;

    if (property_ids.find(i) != property_ids.end()) continue;

    if (id_caches->size() && id_caches->find(i) == id_caches->end())
      continue;

    if ((*id_caches)[i]->size() < min_cand) {
      min_cand = (*id_caches)[i]->size();
      min_cand_root = i;
    }
  }
  return min_cand_root;
}

int TopKSearchPlan::SelectRootRandom(shared_ptr<BGPQuery> &bgp_query,
                                     IDCachesSharePtr &id_caches,
                                     set<int> &property_ids) {
  std::random_device rd;
  std::mt19937 g(rd());
  vector<decltype(total_vars_num_)> root_nodes;

  // get the root id by arg min( candidates * tree depth )
// if all nodes don't have candidates, choose the min tree depth
  for(decltype(total_vars_num_)  i = 0; i< total_vars_num_; i++)
  {
    auto var_descriptor = bgp_query->get_vardescrip_by_index(i);
    if(var_descriptor->var_type_ == VarDescriptor::VarType::Predicate)
      continue;

    if(property_ids.find(i)!=property_ids.end()) continue;

    if(id_caches->size()&&id_caches->find(i)==id_caches->end())
      continue;
    root_nodes.push_back(i);
  }
  std::shuffle(root_nodes.begin(), root_nodes.end(), g);
  return root_nodes[0];
}


void TopKSearchPlan::BuildPlanFromRoot(int min_score_root) {// construct the query tree
  bool*  vars_used_vec = new bool[total_vars_num_];
  memset(vars_used_vec, 0, total_vars_num_ *sizeof(bool));
  auto root_id = min_score_root;
#ifdef SetRootVar0
  root_id = 0;
#endif
  vars_used_vec[root_id] = true;
  stack<TopKTreeNode*> explore_id;

  TopKTreeNode* r = new TopKTreeNode;
  r->var_id = root_id;
  id_node_mapping_[r->var_id] = r;
  var_child_order_[r->var_id] = 0;
  explore_id.push(r);

  vector<StepOperation> edges_after_;
  while (!explore_id.empty()) {
    auto now_node = explore_id.top();
    auto now_id = now_node->var_id;
    explore_id.pop();
    unsigned child_num = neighbours_[now_id].size();
    for (unsigned child_i = 0; child_i < child_num; ++child_i) {
      auto child_id = neighbours_[now_id][child_i];
      // below are vectors
      auto two_var_predicate_ids = predicates_ids_[now_id][child_i];
      auto two_var_predicate_constants = predicates_constant_[now_id][child_i];
      auto two_var_directions = directions_[now_id][child_i];

      if(vars_used_vec[child_id])
        continue;

      vars_used_vec[child_id] = true;
      auto child_tree = new TopKTreeNode;
      child_tree->var_id = child_id;
      now_node->descendents_.push_back(child_tree);
      auto tree_edge_ptr = make_shared<TopKPlanUtil::TreeEdge>();
      tree_edge_ptr->predicate_constant_ = move(two_var_predicate_constants);
      tree_edge_ptr->predicate_ids_ = move(two_var_predicate_ids);
      tree_edge_ptr->directions_ = move(two_var_directions);
      now_node->tree_edges_.push_back(tree_edge_ptr);
      tree_edge_ptr->ChangeOrder();
      explore_id.push(child_tree);
      id_node_mapping_[child_id] = child_tree;
    }
  }
  delete[] vars_used_vec;
  tree_root_ = r;
  AdjustOrder();
  max_leave_distance_.reset(new size_t[total_vars_num_ + 1]);
  BuildMaxDistanceFromLeaf(tree_root_);
}

void TopKSearchPlan::GetPlanCostModel(shared_ptr<BGPQuery> bgp_query, KVstore *kv_store, Statistics *statistics, const QueryTree::Order& expression,
                                      IDCachesSharePtr id_caches,int k,
                                      std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){


//  auto first_cycle = this->FindCycle();
//  if(!first_cycle.empty())
//    this->CutCycle(bgp_query,kv_store,statistics,id_caches);
  set<int> property_ids;
  for(auto& var_co_pair:*var_coefficients)
    if(var_co_pair.second!=0)
      property_ids.insert(bgp_query->get_var_id_by_name(var_co_pair.first));

  int root_id;
#ifdef SetRootVar0
  root_id = 0;
#endif
  switch (global_root_strategy) {
    case RootStrategy::Ours:
      //cerr<<101<<endl;
      root_id = this->SelectRootCostModelOurs(bgp_query, id_caches, total_vars_num_, k,property_ids);
      break;
    case RootStrategy::MinCand:
      //cerr<<102<<endl;
      root_id = this->SelectRootMinCand(bgp_query,id_caches,property_ids);
      break;
    case RootStrategy::MinDepth:
      //cerr<<103<<endl;
      root_id = this->SelectRootMinDepth(bgp_query,id_caches,property_ids);
      break;
    case RootStrategy::Random:
      //cerr<<104<<endl;
      root_id = this->SelectRootRandom(bgp_query,id_caches,property_ids);
      break;
  }
  cout << "choose root["<<root_id<<"] "<<bgp_query->get_var_name_by_id(root_id)<<endl;
  root_var_name = bgp_query->get_var_name_by_id(root_id);
  BuildPlanFromRoot(root_id);
}

TopKSearchPlan::~TopKSearchPlan() {
  if(tree_root_== nullptr)
    return;
  std::stack<TopKTreeNode*> to_delete;
  to_delete.push(this->tree_root_);
  while(!to_delete.empty())
  {
    auto tree_node = to_delete.top();
    to_delete.pop();
    for(const auto child: tree_node->descendents_)
    {
      to_delete.push(child);
    }
    delete tree_node;
  }
}

void TopKSearchPlan::DebugInfo(shared_ptr<BGPQuery> bgp_query, KVstore *kv_store) {

  std::stack<decltype(this->tree_root_)> h;
  h.push(this->tree_root_);
  while(!h.empty())
  {
    auto now_node = h.top();
    h.pop();
    auto child_num = now_node->descendents_.size();
#ifdef TOPK_DEBUG_INFO
    cout << "pop "  <<bgp_query->get_var_name_by_id( now_node->var_id);
    cout<<"["<<now_node->var_id << "] children num:" << child_num<<endl;
#endif
    for(decltype(child_num) child_i=0;child_i<child_num;child_i++)
    {
      auto &child_p = now_node->descendents_[child_i];
      auto &tree_edge_ptr = now_node->tree_edges_[child_i];
      auto child_id = child_p->var_id;

#ifdef TOPK_DEBUG_INFO
      cout << "\t child "<<bgp_query->get_var_name_by_id(child_id)<<"["<<child_id<<"]"<<endl;
      auto edge_num = tree_edge_ptr->predicate_constant_.size();
      for(decltype(edge_num) j = 0;j<edge_num;j++)
      {
        cout << "\t \t edge-------------"<<j<<"----------------"<<endl;
        cout << "\t \t direction:";
        if(tree_edge_ptr->directions_[j]==TopKPlanUtil::EdgeDirection::IN)
          cout<<"IN";
        else
          cout<<"OUT";
        cout <<endl<< "\t \t predicate_constant:"<<tree_edge_ptr->predicate_constant_[j];
        cout <<endl<< "\t \t predicate:";
        if(tree_edge_ptr->predicate_constant_[j])
          cout<<kv_store->getPredicateByID(tree_edge_ptr->predicate_ids_[j]);
        else
          cout<<bgp_query->get_var_name_by_id( tree_edge_ptr->predicate_ids_[j])<<"["
              <<tree_edge_ptr->predicate_ids_[j]<<"]";
        cout<<endl;
      }
      cout << "\t push " << child_p->var_id<<endl;
#endif
      h.push(child_p);
    }

  }
  flush(cout);
  return;
}

/**
 * Change The order so that each node's descendents are OW first and FRs last
 * Making Order Retrieve easier
 */
void TopKSearchPlan::AdjustOrder() {
  // count fr and ow information
  // fr[i] means node i will corresponds to a fr in its
  // parent node
  bool is_fr[this->total_vars_num_+1];
  for(size_t i = 0; i < this->total_vars_num_+1; ++i)
    is_fr[i] = false;

  stack<TopKTreeNode*> node_to_visit;
  node_to_visit.push(this->tree_root_);

  while (!node_to_visit.empty()) {
    auto top = node_to_visit.top();
    top->descendents_fr_num_ =0 ;
    top->descendents_ow_num_ =0 ;
    if(!top->descendents_.empty())
      is_fr[top->var_id] = true;
    node_to_visit.pop();
    for (auto child_p : top->descendents_)
      node_to_visit.push(child_p);
  }
  // now to adjust the order
  node_to_visit.push(this->tree_root_);

  while (!node_to_visit.empty()) {
    auto top = node_to_visit.top();

    // [0,const_end) stores constant edge
    // [const_end,end) stores variable edge
    size_t ow_end = 0;
    unsigned descendents_num = top->descendents_.size();
    for (unsigned i = 0; i < descendents_num; ++i) {
      auto i_node_var = top->descendents_[i]->var_id;
      // move it to the left part
      if (!is_fr[i_node_var]) {
        // exchange i with const_end
        if (i != ow_end) {
          std::swap(top->tree_edges_[i],top->tree_edges_[ow_end]);
          std::swap(top->descendents_[i], top->descendents_[ow_end]);
        }
        ow_end++;
      }
    }

    top->descendents_ow_num_ = ow_end;
    top->descendents_fr_num_ = descendents_num - ow_end;

    node_to_visit.pop();
    for (unsigned int child_counter = 0; child_counter < top->descendents_.size(); ++child_counter) {
      auto child_p = top->descendents_[child_counter];
      auto child_id = child_p->var_id;
      if(is_fr[child_id])
        node_to_visit.push(child_p);
      this->var_child_order_[child_id] = child_counter;
    }
  }
}

/**
 * Using Topological method to find if the query has a cycle
 * @return the found cycle. Empty if not found
 */
std::vector<int> TopKSearchPlan::FindCycle() {
  std::vector<int> result_cycle;
  auto degrees = map<int,size_t>();
  stack<int> one_degree;

  for(auto &pair:this->neighbours_)
  {
    auto var_id = pair.first;
    auto var_degree = pair.second.size();
    if(var_degree>1)
      degrees[var_id] = var_degree;
    if(var_degree==1)
      one_degree.push(var_id);
  }

  while(!one_degree.empty())
  {
    auto top_id = one_degree.top();
    one_degree.pop();
    auto& top_neighbours = this->neighbours_[top_id];
    for(auto neighbour_id:top_neighbours)
    {
      if(degrees.find(neighbour_id)==degrees.end())
        continue;
      if(--degrees[neighbour_id]==1) {
        one_degree.push(neighbour_id);
        degrees.erase(neighbour_id);
      }
    }
  }

  if(degrees.empty())
    return result_cycle;

  // a cycle exists in the remaining 'degrees', find it
  set<int> possible_vars;
  for(auto& pair:degrees)
    possible_vars.insert(pair.first);

  auto cycle_start = *possible_vars.begin();
  set<int> walk_pass_vars;
  walk_pass_vars.insert(cycle_start);
  result_cycle.push_back(cycle_start);
  bool found = walk(possible_vars,walk_pass_vars,result_cycle);
  if(!found)
    throw string(" Not Expected Situation for topology");
  return move(result_cycle);
}

/**
 * To judge if the query graph is too
 * complex (more than 2 loops ) to
 * use top k algorithm
 * @return if this query graph worth of top-k
 */
bool TopKSearchPlan::SuggestTopK()
{
  auto cycle = this->FindCycle();
  if(cycle.empty())
    return true;
  auto id1 = cycle[0];
  auto id2 = cycle[1];
  auto id1_neighbour_backup = this->neighbours_[id1];
  auto id2_neighbour_backup = this->neighbours_[id2];

  auto &old_id1_neighbour = this->neighbours_[id1];
  auto &old_id2_neighbour = this->neighbours_[id2];
  old_id1_neighbour = vector<TYPE_ENTITY_LITERAL_ID>();
  old_id2_neighbour = vector<TYPE_ENTITY_LITERAL_ID>();

  for(auto id1_nei: id1_neighbour_backup)
    if(id1_nei!=id2)
      old_id1_neighbour.push_back(id1_nei);

  for(auto id2_nei: id2_neighbour_backup)
    if(id2_nei!=id1)
      old_id2_neighbour.push_back(id2_nei);

  auto second_cycle = this->FindCycle();
  bool worth_try = second_cycle.empty();

  this->neighbours_[id1] = id1_neighbour_backup;
  this->neighbours_[id2] = id2_neighbour_backup;
  return worth_try;
}

/**
 * A dfs way to find cycle
 */
bool TopKSearchPlan::walk(set<int> &possible_vars,set<int> &walk_pass_vars, vector<int> &result_cycle)
{
  auto var_now = result_cycle.back();
  auto& next_vars = this->neighbours_[var_now];
  for(auto next_var:next_vars)
  {
    if(possible_vars.find(next_var)==possible_vars.end())
      continue;

    // We have found a cycle
    if(walk_pass_vars.find(next_var)!=walk_pass_vars.end())
    {
      auto start_it = result_cycle.begin();
      while(*start_it!=next_var)
        start_it++;
      vector<int> tem(start_it,result_cycle.end());
      result_cycle = tem;
      return true;
    }

    // if not found, walk deeper
    result_cycle.push_back(next_var);
    walk_pass_vars.insert(next_var);
    bool found = walk(possible_vars,walk_pass_vars,result_cycle);
    if(found)
      return found;
    result_cycle.pop_back();
    walk_pass_vars.erase(next_var);
    // to try next
  }
  return false;
}

/**
 * Cut A edge from the cycle, so that this query
 * will be processed as tree query. After each
 * result, we will check if the edge exist
 * @return if we cut an edge
 */
bool TopKSearchPlan::CutCycle(shared_ptr<BGPQuery> bgp_query, KVstore *kv_store, Statistics *statistics,
                              IDCachesSharePtr id_caches) {
  auto cycle = this->FindCycle();
  if(cycle.empty())
    return false;
  this->is_cycle_graph_ = true;
  vector<double> selectivity(cycle.size());
  auto edge_num = cycle.size();
  cycle.push_back(cycle[0]);
  decltype(edge_num) choose_one = 0;
  // choose a max selectivity
  double max_selectivity = 1.1;
  for(decltype(edge_num) i=0;i<edge_num;i++)
  {
    // estimate selectivity of a->b
    auto a_id = cycle[i];
    auto b_id = cycle[i+1];
    auto a_var = bgp_query->get_vardescrip_by_index(a_id);
    auto b_string = bgp_query->get_var_name_by_id(b_id);
    auto &edges = a_var->so_edge_index_;
    double min_sel = 1;

    shared_ptr<IDList> a_cache;
    shared_ptr<IDList> b_cache;

    auto a_it = id_caches->find(a_id);
    if(a_it==id_caches->end())
      continue;
    auto b_it = id_caches->find(b_id);
    if(b_it==id_caches->end())
      continue;

    a_cache = a_it->second;
    b_cache = b_it->second;

    for(auto edge_id:edges) {
      auto triple = bgp_query->get_triple_by_index(edge_id);
      auto& s_string = triple.subject;
      auto& p_string = triple.predicate;
      auto& o_string = triple.object;
      if(s_string!=b_string && o_string != b_string)
        continue;
      auto predicate_constant = p_string[0] != '?';
      TYPE_PREDICATE_ID predicate_id = -1;
      if(predicate_constant)
        predicate_id = kv_store->getIDByPredicate(p_string);

      double sel;
      if(s_string==b_string) {
        sel = PlanGenerator::estimate_one_edge_selectivity(predicate_id,predicate_constant,kv_store,
                                                     b_cache,a_cache);
      }
      else // s_string==a_string
      {
        sel = PlanGenerator::estimate_one_edge_selectivity(predicate_id,predicate_constant,kv_store,
                                                     a_cache,b_cache);
      }
      min_sel = std::min(min_sel,sel);
    }

    if(min_sel>max_selectivity)
    {
      max_selectivity = min_sel;
      choose_one = i;
    }
  }

  // Now cut the edge c-d
  // and record the edge in the class member
  auto c_id = cycle[choose_one];
  auto d_id = cycle[choose_one+1];
  DeleteEdge(c_id,d_id);
  DeleteEdge(d_id,c_id);
  auto c_var = bgp_query->get_vardescrip_by_index(c_id);
  auto d_string = bgp_query->get_var_name_by_id(d_id);
  auto &edges = c_var->so_edge_index_;

  for(auto edge_id:edges) {
    auto triple = bgp_query->get_triple_by_index(edge_id);
    auto &s_string = triple.subject;
    auto &p_string = triple.predicate;
    auto &o_string = triple.object;
    if (s_string != d_string && o_string != d_string)
      continue;
#ifdef TOPK_DEBUG_INFO
    cout<<"cut edge "<<s_string<<" "<<p_string<<" "<<o_string<<endl;
#endif
    auto predicate_constant = p_string[0] != '?';
    TYPE_PREDICATE_ID predicate_id = -1;
    if (predicate_constant)
      predicate_id = kv_store->getIDByPredicate(p_string);

    EdgeInfo edge_info;
    edge_info.p_ = predicate_id;
    if (s_string == d_string) {
      edge_info.s_ = d_id;
      edge_info.o_ = c_id;
    } else {
      edge_info.s_ = c_id;
      edge_info.o_ = d_id;
    }
    edge_info.join_method_ = JoinMethod::so2p;
    EdgeConstantInfo edge_constant_info(false, predicate_constant, false);
    this->non_tree_edges_.edge_filter_ = make_shared<FeedOneNode>();
    this->non_tree_edges_.edge_filter_->edges_->push_back(edge_info);
    this->non_tree_edges_.edge_filter_->edges_constant_info_->push_back(edge_constant_info);
  }
  this->non_tree_edges_.edge_filter_->node_to_join_ = c_id;
  this->non_tree_edges_.join_type_ = StepOperation::JoinType::EdgeCheck;
  return true;
}


void TopKSearchPlan::DeleteEdge(TYPE_ENTITY_LITERAL_ID a, TYPE_ENTITY_LITERAL_ID b) {
  auto a_nei = neighbours_[a].begin();
  auto a_nei_end = neighbours_[a].end();
  auto a_pre_c = predicates_constant_[a].begin();
  auto a_pre = predicates_ids_[a].begin();
  auto a_direction = directions_[a].begin();
  while(a_nei!=a_nei_end)
  {
    if(*a_nei == b)
    {
      neighbours_[a].erase(a_nei);
      predicates_constant_[a].erase(a_pre_c);
      predicates_ids_[a].erase(a_pre);
      directions_[a].erase(a_direction);
      break;
    }
    a_nei++;
    a_pre_c++;
    a_pre++;
    a_direction++;
  }
}

bool TopKSearchPlan::BuildMaxDistanceFromLeaf(TopKTreeNode *node) {
  auto var_id = node->var_id;
  if (node->descendents_.empty())  // this is leaf
    this->max_leave_distance_[var_id] = 0;
  size_t distance = 0;
  for (auto child : node->descendents_) {
    BuildMaxDistanceFromLeaf(child);
    auto child_distance = this->max_leave_distance_[child->var_id];
    distance = max(distance,child_distance + 1);
  }
  this->max_leave_distance_[var_id] = distance;
  return true;
}

/**
 * Esitimate the cost of build/explore a var's candidate node
 * @param bgp_query
 * @param id_caches
 * @param visit_now
 * @param vars_used_vec
 * @param total_vars_num
 * @return
 */
double TopKSearchPlan::CostEstimate(shared_ptr<BGPQuery> bgp_query,
                                    IDCachesSharePtr id_caches,
                                    TYPE_ENTITY_LITERAL_ID visit_now,
                                    vector<bool> vars_used_vec,
                                    std::size_t total_vars_num,
                                    int k){
  double cost = 0;
  vars_used_vec[visit_now] = true;
  auto &children = this->neighbours_[visit_now];
  if(children.size()>1)
  {
    auto HeapFlowSize = this->GetMinHeapNode(k,children.size());
    cost += HeapFlowSize + k*log2(HeapFlowSize);
  }
  for(const auto& child_id: this->neighbours_[visit_now]) {
    if (!vars_used_vec[child_id]) {
      // for numeric variables, we don't have estimate, so ignore
      // cost += 0
      if(id_caches->find(child_id)==id_caches->end())
        continue;
      auto now_size = id_caches->find(visit_now)->second->size();
      auto descend_cost = this->CostEstimate(bgp_query, id_caches, child_id, vars_used_vec, total_vars_num,k);
      auto descend_size = id_caches->find(child_id)->second->size();
      double ratio = std::max((double)1.0, 1.0*(double)descend_size / now_size);
      cost += descend_cost + log2(ratio)*k;
      vars_used_vec[child_id] = true;
    }
  }
  return cost;
}
int TopKSearchPlan::SelectRootCostModelOurs(shared_ptr<BGPQuery> bgp_query,
                                            IDCachesSharePtr id_caches,
                                            std::size_t total_vars_num,
                                            int k,set<int> &property_ids) {
  double min_cost;
  bool first_used = false;
  int selected_id = -1;
  vector<bool> var_used_vec(total_vars_num,false);
  for(int i = 0;i<total_vars_num;i++)
  {
    // literal variables, ignore
    if(id_caches->find(i)==id_caches->end())
      continue;
    for(int j = 0; j<total_vars_num;j++)
      var_used_vec[j] = false;
    var_used_vec[i]=true;
    auto cost = this->CostEstimate(bgp_query,id_caches,i,var_used_vec,total_vars_num,k);
    auto root_C = id_caches->find(i)->second->size();
    cost += k * log2(root_C);
    if( (!first_used) || cost < min_cost )
    {
      first_used = true;
      min_cost = cost;
      selected_id = i;
    }
  }
  return selected_id;
}

double TopKSearchPlan::GetMinHeapNode(int k, int m) {
  //  inline static size_t E(unsigned int balls,unsigned int box)
  //    return Combinatorial[box + balls - 1][box - 1];
  int accumulator = 0;
  int balls = 1;
  for(;balls<2048;balls++)
  {
    accumulator+= Combinatorial[balls+m-1][m-1];
    if(accumulator>=k)
      break;
  }
  accumulator -= Combinatorial[balls+m-1][m-1];
  balls--;
  auto extra_states = k - accumulator;
  // *
  // **
  // *** <- balls

  return extra_states/m + accumulator + Combinatorial[balls+m-1][m-1];
}

/**
 * Reorder the TreeEdge structure, making
 * edges with constant edges pop first
 */
void TopKPlanUtil::TreeEdge::ChangeOrder() {
  // [0,const_end) stores constant edge
  // [const_end,end) stores variable edge
  size_t const_end = 0;
  auto edges_num = predicate_constant_.size();
  for (size_t i = 0; i < edges_num; i++) {
    // const
    if (predicate_constant_[i]) {
      // exchange i with const_end
      if (i != const_end) {
        std::swap(predicate_constant_[i],predicate_constant_[const_end]);
        std::swap(predicate_ids_[i],predicate_ids_[const_end]);
        std::swap(directions_[i],directions_[const_end]);
      }
      const_end++;
    }
  }
}