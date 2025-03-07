/*=============================================================================
# Filename: Executor.h
# Author: Yuqi Zhou
# Mail: zhouyuqi@pku.edu.cn
# Last Modified:  2021/8/16.
=============================================================================*/

#include "Executor.h"
using namespace std;

/**
 * Join A node to the exist table
 * @param old_table
 * @param id_caches
 * @param feed_plan
 * @return the new created table
 */
tuple<bool, IntermediateResult> Executor::JoinANode(IntermediateResult old_table,
                                                    const IDCachesSharePtr id_caches,
                                                    const shared_ptr<FeedOneNode> feed_plan)
{
  auto new_id = feed_plan->node_to_join_;
  auto table_content_ptr = old_table.values_;
  auto new_table = IntermediateResult::OnlyPositionCopy(old_table);
  auto new_records = new_table.values_;
  auto id_pos_map = new_table.id_pos_map;
  auto pos_id_map = new_table.pos_id_map;
  auto edges_info = feed_plan->edges_;
  size_t new_pos = id_pos_map->size();
  (*id_pos_map)[new_id] = new_pos;
  (*pos_id_map)[new_pos] = new_id;

  /* : each record */
  for (auto record_iterator = table_content_ptr->begin();
       record_iterator != table_content_ptr->end();
       record_iterator++) {
    shared_ptr<IDList> record_candidate_list = ExtendRecordOneNode(feed_plan,
                                                                   id_pos_map,
                                                                   id_caches,
                                                                   new_id,
                                                                   record_iterator);

    /* write to the new table */
    auto record = *record_iterator;
    for (auto new_element:*(record_candidate_list->getList())) {
      auto new_record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>(*record);
      new_record->push_back(new_element);
      new_records->push_back(std::move(new_record));
    }
  }
  return make_tuple(true, new_table);
}

/**
 *
 * @param feed_plan
 * @param is_entity
 * @param is_literal
 * @param is_predicate
 * @param id_caches
 * @return
 */
std::tuple<bool,IntermediateResult> Executor::InitialTableOneNode(const std::shared_ptr<FeedOneNode> feed_plan,
                                                                  bool is_entity,
                                                                  bool is_literal,
                                                                  bool is_predicate,
                                                                  const IDCachesSharePtr id_caches)
{
  IntermediateResult init_table;
  auto &records = init_table.values_;
  auto node_added = feed_plan->node_to_join_;
  init_table.AddNewNode(node_added);
  auto id_position_map = init_table.id_pos_map;
  std::shared_ptr<IDList> first_candidate_list = nullptr;

  if(feed_plan->edges_->size()!=0) {
    // we assert that the plan will not visit contents from empty record
    auto empty_record = records->begin();
    first_candidate_list = ExtendRecordOneNode(feed_plan, id_position_map, id_caches,
                                               node_added, empty_record);
    records = ConvertToTable(first_candidate_list);
    return make_tuple(true, init_table);
  }

  auto id_cache_it = id_caches->find(node_added);
  if (id_cache_it != id_caches->end()) {
    records = ConvertToTable(id_cache_it->second);
    return make_tuple(true, init_table);
  }
  else // No Constant Constraint, Return All IDs
  {
    if(is_predicate)
      records = get<1>(this->GetAllPreId());
    else if(is_entity)
    {
      records = get<1>(this->GetAllSubObjId(is_literal));
    }
    return make_tuple(true, init_table);
  }
}


std::tuple<bool,IntermediateResult> Executor::InitialTableTwoNode(const std::shared_ptr<FeedTwoNode> join_plan,
                                                                  const IDCachesSharePtr id_caches)
{
  IntermediateResult init_table;

  auto records = init_table.values_;
  auto id1 = join_plan->node_to_join_1_;
  auto id2 = join_plan->node_to_join_2_;
  init_table.AddNewNode(id1);
  init_table.AddNewNode(id2);
  auto id_position_map = init_table.id_pos_map;

  // we assert that the plan will not visit contents from empty record
  auto empty_record = records->begin();
  auto record_two_node_list = ExtendRecordTwoNode(join_plan,
                                                  id_position_map,
                                                  id_caches,
                                                  id1,
                                                  id2,
                                                  empty_record);
  auto valid_size = record_two_node_list->size()/2;
  /* write to the new table */
  for (decltype(valid_size) i =0;i<valid_size;i++) {
    auto new_record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    auto id1_this_record = (*record_two_node_list)[i*2];
    auto id2_this_record = (*record_two_node_list)[i*2+1];
    new_record->push_back(id1_this_record);
    new_record->push_back(id2_this_record);
    records->push_back(std::move(new_record));
  }
  return make_tuple(true, init_table);
}



std::tuple<bool, IntermediateResult> Executor::JoinTwoNode(const shared_ptr<FeedTwoNode> join_two_node_,
                                                           IntermediateResult old_table,
                                                          const IDCachesSharePtr id_caches)
{
  auto new_table=IntermediateResult::OnlyPositionCopy(old_table);
  auto id1 = join_two_node_->node_to_join_1_;
  auto id2 = join_two_node_->node_to_join_2_;
  new_table.AddNewNode(id1);
  new_table.AddNewNode(id2);

  auto new_records = new_table.values_;
  auto old_records = old_table.values_;
  auto id_pos_mapping = new_table.id_pos_map;

  /* : each record */
  for (auto record_iterator = old_records->begin();
       record_iterator != old_records->end();
       record_iterator++) {
    auto record_two_node_list = ExtendRecordTwoNode(join_two_node_,
                                                                   id_pos_mapping,
                                                                   id_caches,
                                                                   id1,
                                                                   id2,
                                                                   record_iterator);
    auto valid_size = record_two_node_list->size()/2;
    /* write to the new table */
    auto record = *record_iterator;
    for (decltype(valid_size) i =0;i<valid_size;i++) {
      auto new_record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>(*record);
      auto id1_this_record = (*record_two_node_list)[i*2];
      auto id2_this_record = (*record_two_node_list)[i*2+1];
      new_record->push_back(id1_this_record);
      new_record->push_back(id2_this_record);
      new_records->push_back(std::move(new_record));
    }
  }

  return make_tuple(true, new_table);
}


/**
 * join two table with vars in one_step_join_table
 * @param join_plan: define join vars
 * @param table_a
 * @param table_b
 * @return new table, columns are made up of
 * [ big table vars ][ small table vars - common vars]
 */
tuple<bool,IntermediateResult> Executor::JoinTable(const shared_ptr<JoinTwoTable> join_plan,
                                                   IntermediateResult table_a,
                                                   IntermediateResult table_b)
{
  auto records_a = table_a.values_;
  auto table_a_id_pos = table_a.id_pos_map;
  auto table_a_pos_id = table_a.pos_id_map;
  auto records_b = table_b.values_;
  auto table_b_id_pos = table_b.id_pos_map;
  auto table_b_pos_id = table_b.pos_id_map;

  IntermediateResult result_table;

  long t1 = Util::get_cur_time();
  if(records_a->empty() || records_b->empty())
    return make_tuple(false,result_table);

  auto new_position_id_mapping = make_shared<PositionValue>();

  auto join_nodes = join_plan->public_variables_;

  // build index in big table
  auto& big_table = records_a->size() > records_b->size() ? records_a : records_b;
  auto& big_id_pos = records_a->size() > records_b->size() ? table_a_id_pos : table_b_id_pos;
  auto& big_pos_id = records_a->size() > records_b->size() ? table_a_pos_id : table_b_pos_id;
  auto& small_table = records_a->size() <= records_b->size() ? records_a : records_b;
  auto& small_id_pos = records_a->size() <= records_b->size() ? table_a_id_pos : table_b_id_pos;
  auto& small_pos_id = records_a->size() <= records_b->size() ? table_a_pos_id : table_b_pos_id;

  auto indexed_result = unordered_map<
      /*key*/ vector<TYPE_ENTITY_LITERAL_ID>,
      /*value*/ shared_ptr<vector<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>,
      /*hash function*/ container_hash<vector<TYPE_ENTITY_LITERAL_ID>>
  >();

  // the new order [ big table vars ][ small table vars - common vars]
  vector<TYPE_ENTITY_LITERAL_ID> common_variables_position_big;
  vector<TYPE_ENTITY_LITERAL_ID> common_variables_position_small;
  for(auto common_variable:*(join_plan->public_variables_))
  {
    common_variables_position_big.push_back((*big_id_pos)[common_variable]);
    common_variables_position_small.push_back((*small_id_pos)[common_variable]);
  }
  auto &public_vars = join_plan->public_variables_;

  auto big_table_record_len = big_id_pos->size();
  for(decltype(big_table_record_len) i =0;i<big_table_record_len;i++)
    (*new_position_id_mapping)[i]= (*big_pos_id)[i];

  auto small_table_record_len = small_pos_id->size();
  for(decltype(small_table_record_len) i =0;i<small_table_record_len;i++) {
    auto small_id = (*small_pos_id)[i];
    if(find(public_vars->begin(),public_vars->end(),small_id)!=public_vars->end())
      continue;
    (*new_position_id_mapping)[new_position_id_mapping->size()] = small_id;
  }

  auto common_variables_size = join_plan->public_variables_->size();
  for(const auto& big_record:*(big_table))
  {
    vector<TYPE_ENTITY_LITERAL_ID> result_index(common_variables_size);
    for(auto common_position:common_variables_position_big)
    {
      result_index.push_back((*big_record)[common_position]);
    }
    if(indexed_result.find(result_index)!=indexed_result.end())
    {
      indexed_result[result_index]->push_back(big_record);
    }
    else
    {
      auto tmp_vector = make_shared<vector<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();
      tmp_vector->push_back(big_record);
      indexed_result[result_index] = tmp_vector;
    }
  }


  /* Now Do the Matching , first calculate which ids in table small should be added */
  set<TYPE_ENTITY_LITERAL_ID> public_variables;
  for(auto variable_id:*join_nodes)
  {
    public_variables.insert(variable_id);
  }
  vector<TYPE_ENTITY_LITERAL_ID> small_table_inserted_variables_position;
  for(decltype(small_pos_id->size()) i =0;i<small_pos_id->size();i++)
  {
    auto small_node = (*small_pos_id)[i];
    if (public_variables.find(small_node)!=public_variables.end())
      continue;
    small_table_inserted_variables_position.push_back(i);
  }

  auto result_records = result_table.values_;
  /* do the matching */
  for(const auto& small_record:*small_table)
  {
    vector<TYPE_ENTITY_LITERAL_ID> result_index(common_variables_size);
    for(auto common_position:common_variables_position_small)
      result_index.push_back((*small_record)[common_position]);

    /* the index is in the big table */
    if(indexed_result.find(result_index)!=indexed_result.end())
    {
      vector<TYPE_ENTITY_LITERAL_ID> small_record_inserted;
      for(auto small_position:small_table_inserted_variables_position)
        small_record_inserted.push_back((*small_record)[small_position]);
      auto matched_content = indexed_result[result_index];
      for(const auto& matched_big_record:*matched_content)
      {
        auto result_record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>(*matched_big_record);
        for(auto small_inserted_element:small_record_inserted)
          result_record->push_back(small_inserted_element);
        result_records->push_back(result_record);
      }
    }
    /* if not, ignore */
  }
  result_table.pos_id_map = new_position_id_mapping;
  auto id_pos_map = result_table.id_pos_map;
  for(const auto& pos_id_pair:*new_position_id_mapping)
    (*id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  return make_tuple(true, result_table);
}


/**
 *
 * @param check_plan
 * @param old_table
 * @param id_caches
 * @return a new IntermediateResult, changing the return one will not change the old table
 */
tuple<bool,IntermediateResult> Executor::ANodeEdgesConstraintFilter(shared_ptr<FeedOneNode> check_plan,
                                                                    IntermediateResult old_table,
                                                                    IDCachesSharePtr id_caches) {
  auto new_table = IntermediateResult::OnlyPositionCopy(old_table);
  auto table_content = old_table.values_;
  auto id_pos_mapping =old_table.id_pos_map;

  auto edge_info_vec = check_plan->edges_;
  auto edge_constant_info_vec = check_plan->edges_constant_info_;
  auto edge_num = check_plan->edges_->size();
  for(decltype(edge_num) i=0;i<edge_num;i++)
  {
    auto edge_info = (*edge_info_vec)[i];
    auto edge_constant_info = (*edge_constant_info_vec)[i];
    auto step_result = this->OneEdgeConstraintFilter(edge_info, edge_constant_info, table_content, id_pos_mapping, id_caches);
    if(get<0>(step_result))
      table_content = get<1>(step_result);
  }
  new_table.values_ = table_content;
  return make_tuple(true,new_table);
}

bool Executor::CacheConstantCandidates(shared_ptr<FeedOneNode> one_step, IDCachesSharePtr id_caches) {
  auto node_t = one_step->node_to_join_;
  auto edges = one_step->edges_;
  auto edges_constant = one_step->edges_constant_info_;
  for(decltype(edges->size()) i =0 ;i<edges->size();i++) {
    AddConstantCandidates((*edges)[i],node_t,id_caches);
  }
  return true;
}

/**
 * a simple version of Join::update_answer_list
 * @param valid_id_list
 * @param id_list
 * @param id_list_len
 * @param id_list_prepared : valid_id_list is initialed
 */
void
Executor::UpdateIDList(const shared_ptr<IDList>& valid_id_list, unsigned* id_list, unsigned id_list_len,bool id_list_prepared)
{
  if(id_list_prepared)
  {
    valid_id_list->intersectList(id_list, id_list_len);
  }
  else
  {
    valid_id_list->reserve(id_list_len);
    for(int i = 0; i < id_list_len; i++)
      valid_id_list->addID(id_list[i]);
  }
}

/**
 * This Operation will change id_caches. UpdateNode id_caches with
 * information from step_operation
 * @param step_operation the edge information
 * @param id_caches
 * @return bool
 */
bool Executor::UpdateIDCache(std::shared_ptr<StepOperation> step_operation,IDCachesSharePtr id_caches) {
  //
  auto feed_plan = step_operation->edge_filter_;
  auto var_node_id = feed_plan->node_to_join_;
  auto candidate_list_it = id_caches->find(var_node_id);
  auto candidate_generated = this->CandidatesWithConstantEdge( feed_plan->edges_);
  if(candidate_list_it==id_caches->end())
  {
    (*id_caches)[var_node_id] = candidate_generated;
  }
  else
  {
    auto candidate_already = candidate_list_it->second;
    candidate_already->intersectList(candidate_generated->getList()->data(),candidate_generated->size());
  }
  return true;
}

TableContentShardPtr Executor::ConvertToTable(std::shared_ptr<IDList> id_list) {
  auto result = make_shared<TableContent>();
  auto id_candidate_vec = id_list->getList();
  for (auto var_id: *id_candidate_vec) {
    auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    record->push_back(var_id);
    result->push_back(record);
  }
  return result;
}

/**
 *
 * @param edge_info
 * @param edge_table_info
 * @param table_content_ptr
 * @param id_pos_mapping
 * @param id_caches
 * @return the filtered table. Change the returned one will not infect the old one
 */
tuple<bool,TableContentShardPtr> Executor::OneEdgeConstraintFilter(EdgeInfo edge_info,
                                                                   EdgeConstantInfo edge_table_info,
                                                                   TableContentShardPtr table_content_ptr,
                                                                   PositionValueSharedPtr id_pos_mapping,
                                                                   IDCachesSharePtr id_caches) {
  TYPE_ENTITY_LITERAL_ID var_to_filter= edge_info.getVarToFilter();
  auto var_in_record = id_pos_mapping->find(var_to_filter)!=id_pos_mapping->end();
  if(edge_table_info.ConstantToVar(edge_info))
  {
    TYPE_ENTITY_LITERAL_ID *edge_candidate_list;
    TYPE_ENTITY_LITERAL_ID this_edge_list_len;
    switch (edge_info.join_method_) {
      case JoinMethod::s2p: { // Because if we don't add a pair of '{}', the editor will report a error of redefinition
        auto s_var_constant_id = edge_info.s_;
        this->kv_store_->getpreIDlistBysubID(s_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::s2o: {
        auto s_var_constant_id = edge_info.s_;
        this->kv_store_->getobjIDlistBysubID(s_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2s: {
        auto p_var_constant_id = edge_info.p_;
        this->kv_store_->getsubIDlistBypreID(p_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2o: {
        auto p_var_constant_id = edge_info.p_;
        this->kv_store_->getobjIDlistBypreID(p_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2s: {
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getsubIDlistByobjID(o_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2p: {
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getpreIDlistByobjID(o_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::so2p: {
        auto s_var_constant_id = edge_info.s_;
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getpreIDlistBysubIDobjID(s_var_constant_id,
                                                  o_var_constant_id,
                                                  edge_candidate_list,
                                                  this_edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::sp2o: {
        auto s_var_constant_id = edge_info.s_;
        auto p_var_constant_id = edge_info.p_;
        this->kv_store_->getobjIDlistBysubIDpreID(s_var_constant_id,
                                                  p_var_constant_id,
                                                  edge_candidate_list,
                                                  this_edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::po2s: {
        auto p_var_constant_id = edge_info.p_;
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getsubIDlistByobjIDpreID(o_var_constant_id,
                                                  p_var_constant_id,
                                                  edge_candidate_list,
                                                  this_edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
    }
    auto validate_list = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>(edge_candidate_list,edge_candidate_list+this_edge_list_len);
    delete [] edge_candidate_list;
    return FilterAVariableOnIDList(validate_list,var_to_filter,table_content_ptr,id_pos_mapping);
  }


  auto result_table = make_shared<TableContent>();
  /* : each record */
  for (auto record_iterator = table_content_ptr->begin(); record_iterator != table_content_ptr->end();record_iterator++) {
    TYPE_ENTITY_LITERAL_ID *edge_candidate_list;
    TYPE_ENTITY_LITERAL_ID edge_list_len;
    switch (edge_info.join_method_) {
      case JoinMethod::s2p: { // Because if we don't add a pair of '{}', the editor will report a error of redefinition
        auto s_var_position = (*id_pos_mapping)[edge_info.s_];
        auto s_var_id_this_record = (**record_iterator)[s_var_position];
        this->kv_store_->getpreIDlistBysubID(s_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);

        break;
      }
      case JoinMethod::s2o: {
        auto s_var_position = (*id_pos_mapping)[edge_info.s_];
        auto s_var_id_this_record = (**record_iterator)[s_var_position];
        this->kv_store_->getobjIDlistBysubID(s_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2s: {
        auto p_var_position = (*id_pos_mapping)[edge_info.p_];
        auto p_var_id_this_record = (**record_iterator)[p_var_position];
        this->kv_store_->getsubIDlistBypreID(p_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2o: {
        auto p_var_position = (*id_pos_mapping)[edge_info.p_];
        auto p_var_id_this_record = (**record_iterator)[p_var_position];
        this->kv_store_->getobjIDlistBypreID(p_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2s: {
        auto o_var_position = (*id_pos_mapping)[edge_info.o_];
        auto o_var_id_this_record = (**record_iterator)[o_var_position];
        this->kv_store_->getsubIDlistByobjID(o_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2p: {
        auto o_var_position = (*id_pos_mapping)[edge_info.o_];
        auto o_var_id_this_record = (**record_iterator)[o_var_position];
        this->kv_store_->getpreIDlistByobjID(o_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
        /**
         * Notice: s or o may be constant
         */
      case JoinMethod::so2p: {
        TYPE_ENTITY_LITERAL_ID  s_var_id_this_record;
        TYPE_ENTITY_LITERAL_ID  o_var_id_this_record;

        if(edge_table_info.s_constant_)
        {
          s_var_id_this_record = edge_info.s_;
        }
        else{
          auto s_var_position = (*id_pos_mapping)[edge_info.s_];
          s_var_id_this_record = (**record_iterator)[s_var_position];
        }

        if(edge_table_info.o_constant_)
        {
          o_var_id_this_record = edge_info.o_;
        }
        else{
          auto o_var_position = (*id_pos_mapping)[edge_info.o_];
          o_var_id_this_record = (**record_iterator)[o_var_position];
        }

        this->kv_store_->getpreIDlistBysubIDobjID(s_var_id_this_record,
                                                  o_var_id_this_record,
                                                  edge_candidate_list,
                                                  edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::sp2o: {
        TYPE_ENTITY_LITERAL_ID  s_var_id_this_record;
        TYPE_ENTITY_LITERAL_ID  p_var_id_this_record;
        if(edge_table_info.s_constant_)
        {
          s_var_id_this_record = edge_info.s_;
        }
        else{
          auto s_var_position = (*id_pos_mapping)[edge_info.s_];
          s_var_id_this_record = (**record_iterator)[s_var_position];
        }

        if(edge_table_info.p_constant_)
        {
          p_var_id_this_record = edge_info.p_;
        }
        else{
          auto p_var_position = (*id_pos_mapping)[edge_info.p_];
          p_var_id_this_record = (**record_iterator)[p_var_position];
        }

        this->kv_store_->getobjIDlistBysubIDpreID(s_var_id_this_record,
                                                  p_var_id_this_record,
                                                  edge_candidate_list,
                                                  edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::po2s: {
        TYPE_ENTITY_LITERAL_ID  p_var_id_this_record;
        TYPE_ENTITY_LITERAL_ID  o_var_id_this_record;

        if(edge_table_info.o_constant_)
        {
          o_var_id_this_record = edge_info.o_;
        }
        else{
          auto o_var_position = (*id_pos_mapping)[edge_info.o_];
          o_var_id_this_record = (**record_iterator)[o_var_position];
        }

        if(edge_table_info.p_constant_)
        {
          p_var_id_this_record = edge_info.p_;
        }
        else{
          auto p_var_position = (*id_pos_mapping)[edge_info.p_];
          p_var_id_this_record = (**record_iterator)[p_var_position];
        }

        this->kv_store_->getsubIDlistByobjIDpreID(o_var_id_this_record,
                                                  p_var_id_this_record,
                                                  edge_candidate_list,
                                                  edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      default:
        throw string("OneEdgeConstraintFilter not support join two node");

    }

    /**
     * Now we get the valid candidate list of var_to_filter
     * do a binary to decide whether the record satisfy this edge
     * */
    auto constant_var_id = edge_info.getVarToFilter(edge_table_info);
    bool var_const = constant_var_id.first;
    auto var_id = constant_var_id.second;
    TYPE_ENTITY_LITERAL_ID var_to_filter_id_this_record;
    if(!var_const)
    {
      auto var_to_filter_position = (*id_pos_mapping)[var_to_filter];
      var_to_filter_id_this_record = (**record_iterator)[var_to_filter_position];
    }
    else
      var_to_filter_id_this_record = var_id;

    if(binary_search(edge_candidate_list,edge_candidate_list+edge_list_len,var_to_filter_id_this_record))
    {
      result_table->push_back(*record_iterator);
    }
    delete [] edge_candidate_list;
  }

  return make_tuple(true,result_table);

}

tuple<bool,TableContentShardPtr> Executor::FilterAVariableOnIDList(shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>> candidate_list,
                                                                   TYPE_ENTITY_LITERAL_ID var_id,
                                                                   TableContentShardPtr table_content_ptr,
                                                                   PositionValueSharedPtr id_posing_mapping) {

  auto new_intermediate_table = make_shared<TableContent>();
  auto var_position = (*id_posing_mapping)[var_id];

  /* : each record */
  for (const auto& record_sp : *table_content_ptr){
    auto var_to_filter_id_this_record=(*record_sp)[var_position];
    if(binary_search(candidate_list->begin(),candidate_list->end(),var_to_filter_id_this_record))
    {
      new_intermediate_table->push_back(record_sp);
    }
  }

  return make_tuple(true,new_intermediate_table);

}

bool Executor::AddConstantCandidates(EdgeInfo edge_info,TYPE_ENTITY_LITERAL_ID targetID,IDCachesSharePtr id_caches)
{
  TYPE_ENTITY_LITERAL_ID *edge_candidate_list;
  TYPE_ENTITY_LITERAL_ID this_edge_list_len;
  switch (edge_info.join_method_) {
    case JoinMethod::s2p: { // Because if we don't add a pair of '{}', the editor will report a error of redefinition
      auto s_var_constant_id = edge_info.s_;
      this->kv_store_->getpreIDlistBysubID(s_var_constant_id,
                                           edge_candidate_list,
                                           this_edge_list_len,
                                           true,
                                           this->txn_);
      break;
    }
    case JoinMethod::s2o: {
      auto s_var_constant_id = edge_info.s_;
      this->kv_store_->getobjIDlistBysubID(s_var_constant_id,
                                           edge_candidate_list,
                                           this_edge_list_len,
                                           true,
                                           this->txn_);
      break;
    }
    case JoinMethod::p2s: {
      auto p_var_constant_id = edge_info.p_;
      this->kv_store_->getsubIDlistBypreID(p_var_constant_id,
                                           edge_candidate_list,
                                           this_edge_list_len,
                                           true,
                                           this->txn_);
      break;
    }
    case JoinMethod::p2o: {
      auto p_var_constant_id = edge_info.p_;
      this->kv_store_->getobjIDlistBypreID(p_var_constant_id,
                                           edge_candidate_list,
                                           this_edge_list_len,
                                           true,
                                           this->txn_);
      break;
    }
    case JoinMethod::o2s: {
      auto o_var_constant_id = edge_info.o_;
      this->kv_store_->getsubIDlistByobjID(o_var_constant_id,
                                           edge_candidate_list,
                                           this_edge_list_len,
                                           true,
                                           this->txn_);
      break;
    }
    case JoinMethod::o2p: {
      auto o_var_constant_id = edge_info.o_;
      this->kv_store_->getpreIDlistByobjID(o_var_constant_id,
                                           edge_candidate_list,
                                           this_edge_list_len,
                                           true,
                                           this->txn_);
      break;
    }
    case JoinMethod::so2p: {
      auto s_var_constant_id = edge_info.s_;
      auto o_var_constant_id = edge_info.o_;
      this->kv_store_->getpreIDlistBysubIDobjID(s_var_constant_id,
                                                o_var_constant_id,
                                                edge_candidate_list,
                                                this_edge_list_len,
                                                true,
                                                this->txn_);
      break;
    }
    case JoinMethod::sp2o: {
      auto s_var_constant_id = edge_info.s_;
      auto p_var_constant_id = edge_info.p_;
      this->kv_store_->getobjIDlistBysubIDpreID(s_var_constant_id,
                                                p_var_constant_id,
                                                edge_candidate_list,
                                                this_edge_list_len,
                                                true,
                                                this->txn_);
      break;
    }
    case JoinMethod::po2s: {
      auto p_var_constant_id = edge_info.p_;
      auto o_var_constant_id = edge_info.o_;
      this->kv_store_->getsubIDlistByobjIDpreID(o_var_constant_id,
                                                p_var_constant_id,
                                                edge_candidate_list,
                                                this_edge_list_len,
                                                true,
                                                this->txn_);
      break;
    }
  }
  auto validate_list = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>(edge_candidate_list,edge_candidate_list+this_edge_list_len);
  if(id_caches->find(targetID) == id_caches->end())
  {
    auto new_ids = make_shared<IDList>();
    UpdateIDList(new_ids,edge_candidate_list,this_edge_list_len,false);
    (*id_caches)[targetID] = new_ids;
  }
  else
  {
    auto prepared_IDs = (*id_caches->find(targetID)).second;
    UpdateIDList(prepared_IDs,edge_candidate_list,this_edge_list_len,true);
  }
  delete [] edge_candidate_list;

  return true;
}

tuple<bool, TableContentShardPtr> Executor::GetAllSubObjId(bool need_literal)
{
  set<TYPE_ENTITY_LITERAL_ID> ids;
  for (TYPE_PREDICATE_ID i = 0; i < this->limitID_entity_; ++i) {
    auto entity_str = this->kv_store_->getEntityByID(i);
    if(entity_str!="")
      ids.insert(i);
  }
  if(need_literal) {
    for (TYPE_PREDICATE_ID i = Util::LITERAL_FIRST_ID; i < this->limitID_literal_; ++i) {
      auto entity_str = this->kv_store_->getLiteralByID(i);
      if (entity_str != "")
        ids.insert(i);
    }
  }
  auto result = make_shared<TableContent>();
  for(auto var_id: ids)
  {
    auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    record->push_back(var_id);
    result->push_back(record);
  }
  return make_tuple(true,result);
}


std::tuple<bool, TableContentShardPtr> Executor::GetAllPreId()
{
  set<TYPE_ENTITY_LITERAL_ID> ids;
  for (TYPE_PREDICATE_ID i = 0; i < this->limitID_predicate_; ++i) {
    auto pre_str = this->kv_store_->getEntityByID(i);
    if(pre_str!="")
      ids.insert(i);
  }
  auto result = make_shared<TableContent>();
  for(auto var_id: ids)
  {
    auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    record->push_back(var_id);
    result->push_back(record);
  }
  return make_tuple(true,result);
}

std::shared_ptr<IDList> Executor::ExtendRecordOneNode(shared_ptr<FeedOneNode> one_step_join_node_,
                                                      PositionValueSharedPtr id_pos_mapping,
                                                      IDCachesSharePtr id_caches,
                                                      TYPE_ENTITY_LITERAL_ID new_id,
                                                      list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>::iterator record_iterator) const{
  auto record_candidate_list= make_shared<IDList>();
  auto record_candidate_prepared = false;
  auto edges_info = one_step_join_node_->edges_;
  auto edge_constant_info_vec = one_step_join_node_->edges_constant_info_;
  /* for each edge */
  for (int edge_i = 0;edge_i<edges_info->size();edge_i++) {
    auto edge_info = (*edges_info)[edge_i];
    auto edge_constant_info = (*edge_constant_info_vec)[edge_i];
    if (record_candidate_prepared && record_candidate_list->empty())
      continue;
    TYPE_ENTITY_LITERAL_ID *edge_candidate_list;
    TYPE_ENTITY_LITERAL_ID edge_list_len;
/* through this edge get candidate */
    switch (edge_info.join_method_) {
      case JoinMethod::s2p: { // Because if we don't add a pair of '{}', the editor will report a error of redefinition
        TYPE_ENTITY_LITERAL_ID s_var_id_this_record;
        if(edge_constant_info.s_constant_)
          s_var_id_this_record = edge_info.s_;
        else {
          auto s_var_position = (*id_pos_mapping)[edge_info.s_];
          s_var_id_this_record = (**record_iterator)[s_var_position];
        }
        this->kv_store_->getpreIDlistBysubID(s_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::s2o: {
        TYPE_ENTITY_LITERAL_ID s_var_id_this_record;
        if(edge_constant_info.s_constant_)
          s_var_id_this_record = edge_info.s_;
        else {
          auto s_var_position = (*id_pos_mapping)[edge_info.s_];
          s_var_id_this_record = (**record_iterator)[s_var_position];
        }
        this->kv_store_->getobjIDlistBysubID(s_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2s: {
        TYPE_ENTITY_LITERAL_ID p_var_id_this_record;
        if(edge_constant_info.p_constant_)
          p_var_id_this_record = edge_info.p_;
        else {
          auto p_var_position = (*id_pos_mapping)[edge_info.p_];
          p_var_id_this_record = (**record_iterator)[p_var_position];
        }
        this->kv_store_->getsubIDlistBypreID(p_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2o: {
        TYPE_ENTITY_LITERAL_ID p_var_id_this_record;
        if(edge_constant_info.p_constant_)
          p_var_id_this_record = edge_info.p_;
        else {
          auto p_var_position = (*id_pos_mapping)[edge_info.p_];
          p_var_id_this_record = (**record_iterator)[p_var_position];
        }
        this->kv_store_->getobjIDlistBypreID(p_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2s: {
        TYPE_ENTITY_LITERAL_ID o_var_id_this_record;
        if(edge_constant_info.o_constant_)
          o_var_id_this_record = edge_info.o_;
        else {
          auto o_var_position = (*id_pos_mapping)[edge_info.o_];
          o_var_id_this_record = (**record_iterator)[o_var_position];
        }
        this->kv_store_->getsubIDlistByobjID(o_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2p: {
        TYPE_ENTITY_LITERAL_ID o_var_id_this_record;
        if(edge_constant_info.o_constant_)
          o_var_id_this_record = edge_info.o_;
        else {
          auto o_var_position = (*id_pos_mapping)[edge_info.o_];
          o_var_id_this_record = (**record_iterator)[o_var_position];
        }
        this->kv_store_->getpreIDlistByobjID(o_var_id_this_record,
                                             edge_candidate_list,
                                             edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::so2p: {
        TYPE_ENTITY_LITERAL_ID s_var_id_this_record;
        if(edge_constant_info.s_constant_)
          s_var_id_this_record = edge_info.s_;
        else {
          auto s_var_position = (*id_pos_mapping)[edge_info.s_];
          s_var_id_this_record = (**record_iterator)[s_var_position];
        }

        TYPE_ENTITY_LITERAL_ID o_var_id_this_record;
        if(edge_constant_info.o_constant_)
          o_var_id_this_record = edge_info.o_;
        else {
          auto o_var_position = (*id_pos_mapping)[edge_info.o_];
          o_var_id_this_record = (**record_iterator)[o_var_position];
        }

        this->kv_store_->getpreIDlistBysubIDobjID(s_var_id_this_record,
                                                  o_var_id_this_record,
                                                  edge_candidate_list,
                                                  edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::sp2o: {
        TYPE_ENTITY_LITERAL_ID s_var_id_this_record;
        if(edge_constant_info.s_constant_)
          s_var_id_this_record = edge_info.s_;
        else {
          auto s_var_position = (*id_pos_mapping)[edge_info.s_];
          s_var_id_this_record = (**record_iterator)[s_var_position];
        }

        TYPE_ENTITY_LITERAL_ID p_var_id_this_record;
        if(edge_constant_info.p_constant_)
          p_var_id_this_record = edge_info.p_;
        else {
          auto p_var_position = (*id_pos_mapping)[edge_info.p_];
          p_var_id_this_record = (**record_iterator)[p_var_position];
        }
        this->kv_store_->getobjIDlistBysubIDpreID(s_var_id_this_record,
                                                  p_var_id_this_record,
                                                  edge_candidate_list,
                                                  edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::po2s: {
        TYPE_ENTITY_LITERAL_ID o_var_id_this_record;
        if(edge_constant_info.o_constant_)
          o_var_id_this_record = edge_info.o_;
        else {
          auto o_var_position = (*id_pos_mapping)[edge_info.o_];
          o_var_id_this_record = (**record_iterator)[o_var_position];
        }

        TYPE_ENTITY_LITERAL_ID p_var_id_this_record;
        if(edge_constant_info.p_constant_)
          p_var_id_this_record = edge_info.p_;
        else {
          auto p_var_position = (*id_pos_mapping)[edge_info.p_];
          p_var_id_this_record = (**record_iterator)[p_var_position];
        }

        this->kv_store_->getsubIDlistByobjIDpreID(o_var_id_this_record,
                                                  p_var_id_this_record,
                                                  edge_candidate_list,
                                                  edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }

      default: // s2po p2so o2sp
        throw "ExtendRecordOneNode only support add 1 node each time";
    }

    UpdateIDList(record_candidate_list,
                 edge_candidate_list,
                 edge_list_len,
                 record_candidate_prepared);
    delete [] edge_candidate_list;
/* If candidate is not prepared, then we ues the first edge as candidate */
    record_candidate_prepared = true;
  }

  if(id_caches->find(new_id)!=id_caches->end())
  {
    auto caches_ptr = (*(id_caches->find(new_id))).second;
    record_candidate_list->intersectList(caches_ptr->getList()->data(),caches_ptr->size());
  }

  return record_candidate_list;
}

/**
 *
 * @param one_step_join_node_ only 1 edge
 * @param id_pos_mapping
 * @param id_caches
 * @param new_id1
 * @param new_id2
 * @param record_iterator if all edges are constant and need not to access vars in record_iterator
 * record_iterator can be nullptr
 * @return a vector of [id1 id2]_1 [id1 id2]_2  [id1 id2]_3 [id1 id2]_4
 */
std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>>
Executor::ExtendRecordTwoNode(shared_ptr<FeedTwoNode> one_step_join_node_,
                              PositionValueSharedPtr id_pos_mapping,
                              IDCachesSharePtr id_caches,
                              TYPE_ENTITY_LITERAL_ID new_id1,
                              TYPE_ENTITY_LITERAL_ID new_id2,
                              list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>::iterator record_iterator) const{

  auto two_node_list = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();

  auto edge_info = one_step_join_node_->edges_;
  auto edge_constant_info = one_step_join_node_->edges_constant_info_;

  TYPE_ENTITY_LITERAL_ID *edge_candidate_list;
  TYPE_ENTITY_LITERAL_ID edge_list_len;
/* through this edge get candidate */
  switch (edge_info.join_method_) {
    case JoinMethod::s2po: {
      /* if s is constant, it is edge constraint check, and should not appear in the function*/
      TYPE_ENTITY_LITERAL_ID s_var_id_this_record;
      if(edge_constant_info.s_constant_)
        s_var_id_this_record = edge_info.s_;
      else
      {
        auto s_var_position = (*id_pos_mapping)[edge_info.s_];
        s_var_id_this_record = (**record_iterator)[s_var_position];
      }
      this->kv_store_->getpreIDobjIDlistBysubID(s_var_id_this_record,
                                                edge_candidate_list,
                                                edge_list_len,
                                                true,
                                                this->txn_);
      break;
    }
    case JoinMethod::p2so: {
      TYPE_ENTITY_LITERAL_ID p_var_id_this_record;
      if(edge_constant_info.p_constant_)
        p_var_id_this_record = edge_info.p_;
      else
      {
        auto p_var_position = (*id_pos_mapping)[edge_info.p_];
        p_var_id_this_record = (**record_iterator)[p_var_position];
      }
      this->kv_store_->getsubIDobjIDlistBypreID(p_var_id_this_record,
                                                edge_candidate_list,
                                                edge_list_len,
                                                true,
                                                this->txn_);
      break;
    }
    case JoinMethod::o2ps: {
      TYPE_ENTITY_LITERAL_ID o_var_id_this_record;
      if(edge_constant_info.o_constant_)
        o_var_id_this_record = edge_info.o_;
      else
      {
        auto o_var_position = (*id_pos_mapping)[edge_info.o_];
        o_var_id_this_record = (**record_iterator)[o_var_position];
      }
      this->kv_store_->getpreIDsubIDlistByobjID(o_var_id_this_record,
                                                edge_candidate_list,
                                                edge_list_len,
                                                true,
                                                this->txn_);
      break;
    }
    default:
      throw "ExtendRecordTwoNode only support add 2 node each time";
  }

  bool id1_cached = false;
  decltype((*id_caches->begin()).second->getList()) id1_caches_ptr;

  bool id2_cached = false;
  decltype(id1_caches_ptr) id2_caches_ptr;

  if(id_caches->find(new_id1)!=id_caches->end())
  {
    id1_cached = true;
    id1_caches_ptr = (*(id_caches->find(new_id1))).second->getList();
  }

  if(id_caches->find(new_id2)!=id_caches->end())
  {
    id2_cached= true;
    id2_caches_ptr = (*(id_caches->find(new_id2))).second->getList();
  }

  edge_list_len = edge_list_len/2;
  for(size_t i =0;i<edge_list_len;i++)
  {
    auto id1 = edge_candidate_list[2*i];
    auto id2 = edge_candidate_list[2*i+1];

    if(id1_cached)
      if(!binary_search(id1_caches_ptr->begin(),id1_caches_ptr->end(),id1))
        continue;

    if(id2_cached)
      if(!binary_search(id2_caches_ptr->begin(),id2_caches_ptr->end(),id2))
        continue;

    two_node_list->push_back(id1);
    two_node_list->push_back(id2);
  }
  delete [] edge_candidate_list;
  return two_node_list;
}

/**
 * use edge info to generate IDList
 * @param edge_info_vector e.g. [sp2o 1 2 _][o2s _ _ 5]
 * @return the IDList produced
 */
shared_ptr<IDList> Executor::CandidatesWithConstantEdge(shared_ptr<vector<EdgeInfo>> edge_info_vector) const {
  auto id_candidate = make_shared<IDList>();
  for(int i = 0; i<edge_info_vector->size(); i++)
  {
    auto edge_info = (*edge_info_vector)[i];
    TYPE_ENTITY_LITERAL_ID *edge_candidate_list;
    TYPE_ENTITY_LITERAL_ID this_edge_list_len;
    cout<<"edge["<<i<<"] \n\t"<< edge_info.toString() << "\n\t join method:"<<JoinMethodToString(edge_info.join_method_)<<endl;
    switch (edge_info.join_method_) {
      case JoinMethod::s2p: { // Because if we don't add a pair of '{}', the editor will report a error of redefinition
        auto s_var_constant_id = edge_info.s_;
        this->kv_store_->getpreIDlistBysubID(s_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::s2o: {
        auto s_var_constant_id = edge_info.s_;
        this->kv_store_->getobjIDlistBysubID(s_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2s: {
        auto p_var_constant_id = edge_info.p_;
        this->kv_store_->getsubIDlistBypreID(p_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::p2o: {
        auto p_var_constant_id = edge_info.p_;
        this->kv_store_->getobjIDlistBypreID(p_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2s: {
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getsubIDlistByobjID(o_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::o2p: {
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getpreIDlistByobjID(o_var_constant_id,
                                             edge_candidate_list,
                                             this_edge_list_len,
                                             true,
                                             this->txn_);
        break;
      }
      case JoinMethod::so2p: {
        auto s_var_constant_id = edge_info.s_;
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getpreIDlistBysubIDobjID(s_var_constant_id,
                                                  o_var_constant_id,
                                                  edge_candidate_list,
                                                  this_edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::sp2o: {
        auto s_var_constant_id = edge_info.s_;
        auto p_var_constant_id = edge_info.p_;
        this->kv_store_->getobjIDlistBysubIDpreID(s_var_constant_id,
                                                  p_var_constant_id,
                                                  edge_candidate_list,
                                                  this_edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
      case JoinMethod::po2s: {
        auto p_var_constant_id = edge_info.p_;
        auto o_var_constant_id = edge_info.o_;
        this->kv_store_->getsubIDlistByobjIDpreID(o_var_constant_id,
                                                  p_var_constant_id,
                                                  edge_candidate_list,
                                                  this_edge_list_len,
                                                  true,
                                                  this->txn_);
        break;
      }
    }
    cout<<"get "<<this_edge_list_len<<" result in this edge "<<endl;
    UpdateIDList(id_candidate,edge_candidate_list,this_edge_list_len,i > 0);
  }
  return id_candidate;
}

