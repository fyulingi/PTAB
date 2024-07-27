
#include "DPPTopKUtil.h"

/**
 * Build Iterators for top k search.
 * @return the final DPPFRIterator
 */
DPPFRIterator *DPPUtil::BuildIteratorTree(const shared_ptr<TopKSearchPlan> tree_search_plan, TopKUtil::Env *env,GlobalQueue *global_queue){

  auto root_plan = tree_search_plan->tree_root_;
  auto root_var = root_plan->var_id;

  // from top to down, build the structure
  auto root_candidate_ids = (*env->id_caches)[root_var];

  set<TYPE_ENTITY_LITERAL_ID> root_candidate;
  auto coefficient_it = env->coefficients->find(env->bgp_query->get_var_name_by_id((root_var)));
  bool has_coefficient = coefficient_it != env->coefficients->end();
  double coefficient = has_coefficient? (*coefficient_it).second:0.0;
  // calculating scores
  shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>> node_score = nullptr;

  // build
  for (auto root_id: *root_candidate_ids->getList()) {
    root_candidate.insert(root_id);
  }
#ifdef TOPK_DEBUG_INFO
  cout<<"ROOT has"<< root_candidate.size()<<" ids"<<endl;
#endif
  if(has_coefficient)
    node_score = GetChildNodeScores(coefficient,root_candidate, false, nullptr, nullptr,
                                    nullptr,env);

  env->subtree_has_coefficient_[root_var] = true;

  // each node's descendents are OW first and FRs last
  std::vector<std::map<TYPE_ENTITY_LITERAL_ID ,std::shared_ptr<DPPFRIterator>>> descendents_FRs;

  std::vector<std::map<TYPE_ENTITY_LITERAL_ID, // parent id
      std::pair<std::shared_ptr<DPPOWIterator>, // its OW
      NodeOneChildVarPredicatesPtr>>> // predicate correspond to the OW
  descendents_OWs;

  descendents_FRs.reserve(root_plan->descendents_fr_num_);
  descendents_OWs.reserve(root_plan->descendents_ow_num_);

  std::vector<bool> ow_has_coefficient;
  ow_has_coefficient.reserve(root_plan->descendents_ow_num_);
  std::set<TYPE_ENTITY_LITERAL_ID > deleted_root_ids;
  auto descendents_num = root_plan->descendents_.size();
  unsigned int no_coefficient_fr_subtree_num = 0;
  for (decltype(descendents_num) descendent_i = 0;descendent_i<descendents_num;descendent_i++ )
  {
    auto descendent_plan = root_plan->descendents_[descendent_i];
    auto tree_edges_plan = root_plan->tree_edges_[descendent_i];
    auto is_fr = descendent_plan->descendents_.size() !=0;

    if(is_fr) {
      // 分情况讨论，是FR还是OW
      auto descendent_iterators = GenerateFRs(root_var, descendent_plan->var_id, tree_edges_plan,
                                              root_candidate, deleted_root_ids, descendent_plan, env,global_queue);
      descendents_FRs.push_back(std::move(descendent_iterators));
      if(!env->subtree_has_coefficient_[descendent_plan->var_id])
        no_coefficient_fr_subtree_num++;
    }
    else
    {
      auto descendent_iterators = GenerateOWs(root_var, descendent_plan->var_id, tree_edges_plan,
                                              root_candidate, deleted_root_ids, env);
      descendents_OWs.push_back(std::move(descendent_iterators));
      auto child_coefficient_it  = env->coefficients->find(env->bgp_query->get_var_name_by_id((descendent_plan->var_id)));
      ow_has_coefficient.push_back(child_coefficient_it !=  env->coefficients->end());
    }
  }

  env->ow_coefficient_counts_[root_var] = ow_has_coefficient;
  // assembling all FQ into one FR
  // constructing children' FQs
  auto child_fqs = AssemblingFrOw(root_var,no_coefficient_fr_subtree_num,root_candidate, node_score,env->k,descendents_FRs,descendents_OWs,global_queue);

  auto empty_pre = make_shared<OnePointPredicateVec>();
  auto fr = new DPPFRIterator();
  // constructing parents' FRs
  for(auto root_id_fq:child_fqs) {
    auto root_id = root_id_fq.first;
    auto fq_pointer = root_id_fq.second;
    DPPUtil::RefineTree(root_var,fq_pointer,OrderedListType::FQ,root_plan,env,global_queue);
    fr->Register(root_id,fq_pointer,empty_pre);
  }

#ifdef TOPK_DEBUG_INFO
  cout<<"ROOT has"<< root_candidate.size()<<"ids "<<", FR has "<<child_fqs.size()<<" FQs"<<endl;
#endif
  return fr;
}
// FQ 可能被多个FR访问
void DPPUtil::RefineTree(int var_id,
                const std::shared_ptr<DPPList>& dpp_list,
                OrderedListType type,
                TopKTreeNode *var_tree_node,
                TopKUtil::Env *env,
                GlobalQueue *global_queue)

{
  if(type==OrderedListType::FR)
  {
    // add fr into all its children' 'parent' field
    auto fr = dynamic_pointer_cast<DPPFRIterator>(dpp_list);
    auto raw_fr = fr.get();
    for(const auto &fq_id_pointer_pair:fr->GetFqMap())
    {
      auto fq_pointer = fq_id_pointer_pair.second;
      fq_pointer->Parents().insert(raw_fr);
      RefineTree(var_id,fq_pointer,OrderedListType::FQ,var_tree_node,env,global_queue);
    }
  }
  else if(type==OrderedListType::FQ)
  {
    auto fq = dynamic_pointer_cast<DPPFQIterator>(dpp_list);
    if(fq->visited_)
      return;
    fq->visited_ = true;
    auto raw_fq_iterator = fq.get();
    auto &ow_coefficient = env->ow_coefficient_counts_[var_id];
    auto ow_size = ow_coefficient.size();
    auto total_size = var_tree_node->descendents_.size();
    for(decltype(ow_size) i = 0;i<ow_size;i++)
    {
      auto ow_iterator =fq->GetChild(i);
      if(ow_coefficient[i])
      {
        auto queue_element = unique_ptr<QueueElement>(new QueueElement);
        queue_element->type_ = QueueElement::Type::OW;
        queue_element->fq_iterator_ = (DPPList *) raw_fq_iterator;
        queue_element->self_dpp_list_ = (DPPList *) (ow_iterator.get());
        queue_element->i_th_type = i;
        queue_element->i_element_.reset(new element);
        queue_element->i_element_->cost = ow_iterator->pool_[0].cost;
        queue_element->i_element_->index = 0;
        queue_element->i_element_->node = ow_iterator->pool_[0].node;
        global_queue->Insert(move(queue_element));
      }
      else
        fq->IncreaseBranchNum();
    }
    for(decltype(ow_size) i = ow_size;i<total_size;i++)
    {
      auto fr_plan = var_tree_node->descendents_[i];
      auto fr_pointer = fq->GetChild(i);
      RefineTree(fr_plan->var_id,fr_pointer,OrderedListType::FR,fr_plan,env,global_queue);
    }
  }
  else{
    // for ow, doing nothing, because we have do this
    // in the FQ part
  }
}



/**
 * Assembling FQ Iterators by OW/FR mappings
 * OW first and FR last
 * @return map<FQ id,pointer to DPPFQIterator>
 */
std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFQIterator>>
DPPUtil::AssemblingFrOw(int fq_var,
                        unsigned int no_coefficient_fr_subtree_num,
                        set<TYPE_ENTITY_LITERAL_ID> &fq_ids,
                        std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>> node_scores, int k,
                        vector<std::map<TYPE_ENTITY_LITERAL_ID, shared_ptr<DPPFRIterator> >> &descendents_FRs,
                        std::vector<
                            std::map<TYPE_ENTITY_LITERAL_ID, // parent id
                                     std::pair<shared_ptr<DPPOWIterator>, // its OW
                                               NodeOneChildVarPredicatesPtr>>> // predicate correspond to the OW
                        &descendents_OWs,
                        GlobalQueue *global_queue) {

  auto ow_size = descendents_OWs.size();
  auto fr_size = descendents_FRs.size();
  auto child_type_num =  ow_size+fr_size;
  std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFQIterator>> id_fqs;
  for(auto fq_id:fq_ids)
  {
    double score=0.0;
    if(node_scores!= nullptr)
      score = (*node_scores)[fq_id];
    auto fq = make_shared<DPPFQIterator>(k, fq_id, child_type_num, score,fq_var);
    for(unsigned i=0;i<no_coefficient_fr_subtree_num;i++)
      fq->IncreaseBranchNum();
    auto raw_fq_iterator = fq.get();
    // assembling fq
    for(decltype(child_type_num) i =0;i<ow_size;i++)
    {
      auto ow_iterator = descendents_OWs[i][fq_id].first;
      fq->Insert(ow_iterator);
      fq->AddOneTypePredicate(descendents_OWs[i][fq_id].second);
    }
    for(decltype(child_type_num) i =0;i<fr_size;i++) {
      auto child_fr = descendents_FRs[i][fq_id];
      fq->Insert(descendents_FRs[i][fq_id]);
      child_fr->SetParent(raw_fq_iterator);
    }
    id_fqs[fq_id] = fq;
    //fq->TryGetNext(k);
#ifdef SHOW_SCORE
    cout<<"FQ["<<fq_id<<"], the min score are "<<fq->pool_[0].cost<<endl;
#endif
  }
  return id_fqs;
}

/**
 * Generate each parents' ids corresponding OW iterators,
 * according to DP-P, Global Queue will be inserted once a OW is created
 * @param parent_var
 * @param child_var
 * @param tree_edges_
 * @param parent_var_candidates
 * @param deleted_parents
 * @param child_tree_node
 * @param env
 * @return
 */
std::map<TYPE_ENTITY_LITERAL_ID, // parent id
std::pair<std::shared_ptr<DPPOWIterator>, // its OW
NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
DPPUtil::GenerateOWs(int parent_var, int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge> tree_edges_,
                     std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                     std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                     TopKUtil::Env *env)
{
  auto coefficient_it = env->coefficients->find(env->bgp_query->get_var_name_by_id((child_var)));
  bool has_coefficient = coefficient_it != env->coefficients->end();
  double coefficient = has_coefficient ? (*coefficient_it).second : 0.0;
  if(has_coefficient)
    env->subtree_has_coefficient_[child_var] = true;
  // record the mapping between parent and children
  // IDListWithAppending not only records the children
  // but also records the edges between them
  std::map<TYPE_ENTITY_LITERAL_ID, shared_ptr<IDListWithAppending>  > parent_child;
  std::map<TYPE_ENTITY_LITERAL_ID, std::set<TYPE_ENTITY_LITERAL_ID> > child_parent;

  std::set<TYPE_ENTITY_LITERAL_ID> child_candidates = TopKUtil::ExtendTreeEdge(parent_var_candidates, child_var,
                                                                               deleted_parents, parent_child,
                                                                               child_parent, tree_edges_, env);

  // calculate children's scores
  std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>> node_score =nullptr;
  if(has_coefficient)
    node_score = GetChildNodeScores(coefficient,child_candidates,true,&deleted_parents,&parent_child,
                                    &child_parent,env);

  for(auto deleted_parent:deleted_parents)
    parent_var_candidates.erase(deleted_parent);
  // Return OW iterators

#ifdef SHOW_SCORE
  cout<<"var["<<child_var<<"] has no child, constructing OW iterators"<<endl;
#endif
  std::vector<TYPE_ENTITY_LITERAL_ID> children_ids;
  std::vector<double> children_scores;

  std::map<TYPE_ENTITY_LITERAL_ID, // parent id
      std::pair<shared_ptr<DPPOWIterator>, // its OW
      NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
                                  result;
  for(auto parent_id:parent_var_candidates)
  {
#ifdef SHOW_SCORE
    cout<<"parent var["<<parent_id<<"] has "<<parent_child[parent_id].size()<<" child, its OW "<<endl;
#endif
    auto ow = make_shared<DPPOWIterator>();
#ifdef TOPK_DEBUG_INFO
    ow->var_id = child_var;
#endif
    auto child_it = parent_child[parent_id]->contents_->cbegin();
    auto child_end =  parent_child[parent_id]->contents_->cend();
    auto ow_predicates = make_shared<NodeOneChildVarPredicates>();
    while(child_it != child_end)
    {
      const auto child_id = child_it->first;
      auto predicate_vec = child_it->second;
#ifdef SHOW_SCORE
      cout<<"["<<child_id<<"]"<<" "<<(*node_score)[child_id];
#endif
      children_ids.push_back(child_id);
      (*ow_predicates)[child_id] = predicate_vec;
      if(has_coefficient)
        children_scores.push_back((*node_score)[child_id]);
      child_it++;
    }
#ifdef SHOW_SCORE
    cout<<endl;
      //cout<<parent_id<<"'s OW ["<<child_var<<"]";
#endif
    if(has_coefficient)
      ow->Insert( env->k,children_ids, children_scores);
    else
      ow->Insert(env->k,children_ids);
    result[parent_id] = make_pair(ow,ow_predicates);
    children_ids.clear();
    children_scores.clear();
  }
  return result;

}


/**
 * extend one edge (and the subtree induced) in the query graph
 * parent - children - descendants
 * This function only deals with a parent - children relation in the query
 * graph. Generate children and the call GenerateFRs to prone the
 * generated candidates ,and assemble them as FQ in the end.
 * @param parent_var
 * @param child_var
 * @param parent_var_candidates
 * @param child_tree_node
 * @param env
 * @return each parent id and their children corresponding FRs/OWs
 */
std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFRIterator>>
DPPUtil::GenerateFRs(int parent_var, int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge> tree_edges_,
                     std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                     std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                     TopKTreeNode *child_tree_node, TopKUtil::Env *env,
                     GlobalQueue *global_queue) {
  auto coefficient_it = env->coefficients->find(env->bgp_query->get_var_name_by_id(child_var));
  bool has_coefficient = coefficient_it != env->coefficients->end();
  double coefficient = has_coefficient ? (*coefficient_it).second : 0.0;
  std::set<TYPE_ENTITY_LITERAL_ID> deleted_parent_ids_this_child;
  env->subtree_has_coefficient_[child_var] = has_coefficient;
  // record the mapping between parent and children
  // IDListWithAppending not only records the children
  // but also records the edges between them
  std::map<TYPE_ENTITY_LITERAL_ID, shared_ptr<IDListWithAppending>  > parent_child;
  std::map<TYPE_ENTITY_LITERAL_ID, std::set<TYPE_ENTITY_LITERAL_ID> > child_parent;

  std::set<TYPE_ENTITY_LITERAL_ID> child_candidates = TopKUtil::ExtendTreeEdge(parent_var_candidates, child_var,
                                                                               deleted_parents, parent_child,
                                                                               child_parent, tree_edges_, env);

  // calculate children's scores
  shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>> node_score = nullptr;
  if(has_coefficient)
    node_score = GetChildNodeScores(coefficient,child_candidates,true,&deleted_parents,&parent_child,
                                    &child_parent,env);

  for(auto deleted_parent:deleted_parents)
    parent_var_candidates.erase(deleted_parent);

  // Return FR iterators
  std::set<TYPE_ENTITY_LITERAL_ID > deleted_children;

  std::vector<bool> ow_has_coefficient_vector;
  ow_has_coefficient_vector.reserve(child_tree_node->descendents_ow_num_);

  std::vector<std::map<TYPE_ENTITY_LITERAL_ID ,std::shared_ptr<DPPFRIterator>>> descendents_FRs;
  descendents_FRs.reserve(child_tree_node->descendents_fr_num_);

  std::vector<
      std::map<TYPE_ENTITY_LITERAL_ID, // parent id
      std::pair<std::shared_ptr<DPPOWIterator>, // its OW
      NodeOneChildVarPredicatesPtr>>> // predicate correspond to the OW
  descendents_OWs;
  descendents_OWs.reserve(child_tree_node->descendents_ow_num_);

  auto descendents_num = child_tree_node->descendents_.size();
  unsigned int no_coefficient_fr_subtree_num = 0;

  for (decltype(descendents_num) descendent_i=0;descendent_i<descendents_num;descendent_i++)
  {
    auto descendent_plan = child_tree_node->descendents_[descendent_i];
    auto descendent_tree_edges_plan = child_tree_node->tree_edges_[descendent_i];
    auto is_fr = descendent_plan->descendents_.size() !=0;
    if(is_fr) {
      auto descendent_iterators = GenerateFRs(child_var,
                                              descendent_plan->var_id,
                                              descendent_tree_edges_plan,
                                              child_candidates,
                                              deleted_children,
                                              descendent_plan, env,
                                              global_queue);
      descendents_FRs.push_back(std::move(descendent_iterators));
      if(env->subtree_has_coefficient_[descendent_plan->var_id])
        env->subtree_has_coefficient_[child_var] =true;
      else
        no_coefficient_fr_subtree_num++;
    }
    else
    {
      auto descendent_iterators = GenerateOWs(child_var,
                                              descendent_plan->var_id,
                                              descendent_tree_edges_plan,
                                              child_candidates,
                                              deleted_children,
                                              env);
      descendents_OWs.push_back(std::move(descendent_iterators));
      auto  descendent_coefficient_it  = env->coefficients->find(env->bgp_query->get_var_name_by_id((descendent_plan->var_id)));
      ow_has_coefficient_vector.push_back(descendent_coefficient_it !=  env->coefficients->end());
      env->subtree_has_coefficient_[child_var] = env->subtree_has_coefficient_[child_var] || env->subtree_has_coefficient_[descendent_plan->var_id];
    }
  }

  env->ow_coefficient_counts_[child_var] = ow_has_coefficient_vector;
  auto should_direct_insert =   ! env->subtree_has_coefficient_[child_var];

  // constructing children' FQs
  auto child_fqs = AssemblingFrOw(child_var,no_coefficient_fr_subtree_num,child_candidates, node_score,env->k,descendents_FRs,descendents_OWs,global_queue);
  if(should_direct_insert)
  {
    for(const auto& id_fq:child_fqs)
    {
      id_fq.second->GetFirst(env->k,global_queue);
    }
  }
  // calculate which parent to be deleted
  for(auto deleted_child_id:deleted_children)
  {
    for(auto parent_id:child_parent[deleted_child_id])
    {
      parent_child[parent_id]->Erase(deleted_child_id);
      if(parent_child[parent_id]->MainKeyNum()==0)
      {
        deleted_parents.insert(parent_id);
        parent_var_candidates.erase(parent_id);
      }
    }
  }

  std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFRIterator>> result;

  // constructing parents' FRs
  // also add predicate information into FRs
  for(auto parent_id:parent_var_candidates) {
    auto fr = make_shared<DPPFRIterator>();
    auto edges_predicates = make_shared<NodeOneChildVarPredicates>();
    for(auto child_id_appending_pair: *(parent_child[parent_id])->contents_)
    {
      auto child_id = child_id_appending_pair.first;
      auto appending_list = child_id_appending_pair.second;
      if(should_direct_insert)
      {
        fr->Insert(child_id,child_fqs[child_id],appending_list);
      }
      else
        fr->Register(child_id,child_fqs[child_id],appending_list);
    }
    result[parent_id] = fr;
    if(should_direct_insert)
      fr->GetFirst();
  }

  return result;

}

void DPPUtil::TriggerSeq(unsigned int k,
                         DPPFQIterator* fq_iterator,
                         DPPList *fr_ow_iterator,
                         unsigned int j_th_type,
                         GlobalQueue *global_queue)
{
  // Different from line 16 in algo.5, maybe a bug in the origin
  // if(fr_ow_iterator->pool_.size()==1)
  if(fq_iterator->EPoolSize()==0)
    fq_iterator->TryTriggeringRootSeq(j_th_type,global_queue);
  else
    fq_iterator->TriggerASingletonSeq(k,j_th_type,global_queue);
}

// line 18-22 in algo.5
void DPPUtil::RorDExpansion(unsigned int k,DPPFQIterator* fq_iterator,
                            GlobalQueue *global_queue,TopKSearchPlan* search_plan)
{
  auto i_th_child = search_plan->GetChildOrder(fq_iterator->var_id);
  if(fq_iterator->pool_.size()==1) {
    for (const auto &t_parent:fq_iterator->Parents()) {
      t_parent->dInsert(k,fq_iterator,i_th_child,global_queue);
    }
  }
  else
  {
    for (const auto &t_parent:fq_iterator->TempParents()) {
      auto index = t_parent.second;
      t_parent.first->rInsert(k,index,fq_iterator,i_th_child,global_queue);
    }
  }
  fq_iterator->CleanTempParents();
}



