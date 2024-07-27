/*=============================================================================
# Filename: Optimizer.cpp
# Author: Yuqi Zhou, Linglin Yang
# Mail: zhouyuqi@pku.edu.cn, fyulingi@gmail.com
=============================================================================*/

#include "Optimizer.h"
#include "../Query/topk/PTAB/DfsUtilDynamic.h"
// #define FEED_PLAN
long UTIL_current_query_time = 0;
vector<double> test_scores{};
vector<vector<unsigned int>> result_ids{};

bool k_value_feed = false;
int benchmark_k = 0;

bool any_k = false;
int any_k_step=0;
vector<long> any_k_time{};

int FQ_NUM=0;
int FR_NUM=0;
int OW_NUM=0;
vector<int> fq_num_vec{};
vector<int> fr_num_vec{};
vector<int> ow_num_vec{};
bool iterator_count_enable=false;
unsigned search_node_count = 0;

Optimizer::Optimizer(KVstore *kv_store,
                     Statistics *statistics,
                     TYPE_TRIPLE_NUM *pre2num,
                     TYPE_TRIPLE_NUM *pre2sub,
                     TYPE_TRIPLE_NUM *pre2obj,
                     TYPE_PREDICATE_ID limitID_predicate,
                     TYPE_ENTITY_LITERAL_ID limitID_literal,
                     TYPE_ENTITY_LITERAL_ID limitID_entity,
                     HopIndex* hop_index,
                     shared_ptr<Transaction> txn
):
    kv_store_(kv_store), statistics(statistics), pre2num_(pre2num),
    pre2sub_(pre2obj),pre2obj_(pre2obj),limitID_predicate_(limitID_predicate),
    limitID_literal_(limitID_literal),limitID_entity_(limitID_entity),hop_index_(hop_index),
    txn_(std::move(txn)), executor_(kv_store,txn,limitID_predicate,limitID_literal,limitID_entity_)
{
/*
  this->current_basic_query_ = -1; // updated by result_list.size()
  this->basic_query_list_= make_shared<vector<shared_ptr<BasicQuery>>>();

  auto basic_query_vector = sparql_query.getBasicQueryVec();
  for(auto basic_query_pointer :basic_query_vector)
  {
    auto basic_query_shared_pointer = make_shared<BasicQuery>(*basic_query_pointer);
    this->basic_query_list_->push_back(basic_query_shared_pointer);
  }

  this->candidate_plans_ = make_shared<vector<tuple<shared_ptr<BasicQuery>, shared_ptr<vector<QueryPlan>>>>>(basic_query_list_->size());

  this->execution_plan_=make_shared<vector<QueryPlan>>();
  this->result_list_=make_shared<vector<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>(); // vector<unsigned*>* result_list;
*/

  // TODO: join_cache_ & cardinality_cache_ not implemented yet.
  this->join_cache_ = make_shared<vector<map<BasicQuery*,vector<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>>>(); // map(sub-structure, result_list)
  this->cardinality_cache_ = make_shared<vector<map<BasicQuery*,shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>>(); // map(sub-structure, cardinality), not in statistics
}

BasicQueryStrategy Optimizer::ChooseStrategy(std::shared_ptr<BGPQuery> bgp_query,QueryInfo *query_info){
  if (!query_info->limit_)
  {
    if(bgp_query->get_triple_num()!=1)
      return BasicQueryStrategy::Normal;
    else if(bgp_query->get_total_var_num() != 3)
      return BasicQueryStrategy::Normal;
    else
      return BasicQueryStrategy::Special;
  }
  else
  {
    if(query_info->ordered_by_expressions_->size()) {
      auto search_plan = make_shared<TopKSearchPlan>(bgp_query, this->kv_store_,this->statistics,
                                                     (*query_info->ordered_by_expressions_)[0],nullptr);
      if(search_plan->SuggestTopK())
        return BasicQueryStrategy::TopK;
      else
        return BasicQueryStrategy::Normal;
    }
    else
      return BasicQueryStrategy::limitK;
  }
}

BasicQueryStrategy Optimizer::ChooseStrategy(BasicQuery *basic_query,QueryInfo *query_info) {
  if (!query_info->limit_) {
    if(basic_query->getTripleNum()!=1)
      return BasicQueryStrategy::Normal;
    else
      return BasicQueryStrategy::Special;
  }
  else
  {
    if(query_info->ordered_by_expressions_->size())
      return BasicQueryStrategy::TopK;
    else
      return BasicQueryStrategy::limitK;
  }
}

std::shared_ptr<std::vector<IntermediateResult>> Optimizer::GenerateResultTemplate(shared_ptr<QueryPlan> query_plan)
{
  auto table_template = make_shared<std::vector<IntermediateResult>>();
  table_template->emplace_back();
  auto &tables = *table_template;
  auto &dfs_operation = *query_plan->join_order_;

  auto join_type = dfs_operation[0].join_type_;
  if(join_type==StepOperation::JoinType::JoinNode)
    tables[0].AddNewNode(dfs_operation[0].join_node_->node_to_join_);
  else if(join_type==StepOperation::JoinType::JoinTwoNodes) {
    tables[0].AddNewNode(dfs_operation[0].join_two_node_->node_to_join_1_);
    tables[0].AddNewNode(dfs_operation[0].join_two_node_->node_to_join_2_);
  }

  auto it = dfs_operation.begin();
  it++;
  while(it!=dfs_operation.end())
  {
    tables.emplace_back( IntermediateResult::OnlyPositionCopy(tables.back()));
    auto join_type = it->join_type_;
    if(join_type==StepOperation::JoinType::JoinNode)
      tables.back().AddNewNode( it->join_node_->node_to_join_);
    else if(join_type==StepOperation::JoinType::JoinTwoNodes) {
      tables.back().AddNewNode( it->join_two_node_->node_to_join_1_);
      tables.back().AddNewNode( it->join_two_node_->node_to_join_2_);
    }
  }
  return table_template;
}

tuple<bool, IntermediateResult> Optimizer::ExecutionDepthFirst(shared_ptr<BGPQuery> bgp_query,
                                                               shared_ptr<QueryPlan> query_plan,
                                                               QueryInfo query_info,
                                                               PositionValueSharedPtr id_pos_mapping) {
  auto limit_num = query_info.limit_num_;
  cout<<"Optimizer::ExecutionDepthFirst query_info.limit_num_="<<query_info.limit_num_<<endl;

  auto var_candidates_cache = make_shared<unordered_map<TYPE_ENTITY_LITERAL_ID,shared_ptr<IDList>>>();
  for(auto& constant_generating_step: *(query_plan->constant_generating_lists_))
  {
    executor_.CacheConstantCandidates(constant_generating_step,var_candidates_cache);
  };
  auto &first_operation = (*query_plan->join_order_)[0];
  tuple<bool, IntermediateResult> step_result;
  IntermediateResult first_table;
  if(first_operation.join_type_==StepOperation::JoinType::JoinNode)
  {
    auto property_array = PrepareInitial(bgp_query,first_operation.join_node_);
    auto is_entity = get<0>(property_array);
    auto is_literal = get<1>(property_array);
    auto is_predicate = get<2>(property_array);
    step_result = executor_.InitialTableOneNode(first_operation.join_node_,is_entity,
                                                is_literal,is_predicate,var_candidates_cache);
    first_table = get<1>(step_result);
  }
  else // Two Nodes
  {
    step_result  = executor_.InitialTableTwoNode(first_operation.join_two_node_,var_candidates_cache);
    first_table = get<1>(step_result);

  }
  if( query_plan->join_order_->size()==1 || first_table.values_->empty())
    return make_tuple(true, first_table);

  auto first_candidates_list = first_table.values_;
  int now_result = 0;
  vector<TableContentShardPtr> result_container;

  auto table_template = GenerateResultTemplate(query_plan);

  while (now_result <= limit_num) {
    auto tmp_result = make_shared<TableContent>();
    tmp_result->push_back(first_candidates_list->front());
    first_candidates_list->pop_front();
    auto first_var_one_point_result =
        this->DepthSearchOneLayer(query_plan, 1, now_result, limit_num, tmp_result,var_candidates_cache,table_template);

    if (get<0>(first_var_one_point_result))
    {
      auto one_point_inter_result = get<1>(first_var_one_point_result).values_;
      if (!one_point_inter_result->empty())
        result_container.push_back(one_point_inter_result);
      if (now_result >= limit_num)
        break;
    }
    if (first_candidates_list->empty())
      break;
  }

  IntermediateResult final_result;
  final_result.pos_id_map = table_template->back().pos_id_map;
  final_result.id_pos_map = table_template->back().id_pos_map;

  cout<<"Optimizer::ExecutionDepthFirst result_container.size()="<<result_container.size()<<endl;
  if(result_container.empty())
    return make_tuple(true,final_result);
  int counter = 0;
  /* merge result */
  for(unsigned int i=1;i<result_container.size();i++)
    for(const auto& record:*(result_container[i])) {
      final_result.values_->push_back(record);
      if( ++counter == limit_num)
        break;
    }
  cout<<"Optimizer::ExecutionDepthFirst final_result.size()="<<final_result.values_->size()<<endl;
  return make_tuple(true,final_result);
}

/**
 *
 * @param query_plan
 * @param layer_count
 * @param result_number_till_now
 * @param limit_number
 * @param table_content_ptr
 * @param id_pos_mapping
 * @param id_caches
 * @return
 */
tuple<bool,IntermediateResult> Optimizer::DepthSearchOneLayer(shared_ptr<QueryPlan> query_plan,
                                                              int layer_count,
                                                              int &result_number_till_now,
                                                              int limit_number,
                                                              TableContentShardPtr table_content_ptr,
                                                              IDCachesSharePtr id_caches,
                                                              shared_ptr<vector<IntermediateResult>> table_template) {
  IntermediateResult old_table;
  old_table.values_ = table_content_ptr;
  old_table.pos_id_map = (*table_template)[layer_count-1].pos_id_map;
  old_table.id_pos_map = (*table_template)[layer_count-1].id_pos_map;

  IntermediateResult layer_result;
  layer_result.pos_id_map = (*table_template)[layer_count].pos_id_map;
  layer_result.id_pos_map = (*table_template)[layer_count].id_pos_map;

  auto one_step = (*(query_plan->join_order_))[layer_count];
  TableContentShardPtr step_table;
  switch (one_step.join_type_) {
    case StepOperation::JoinType::JoinNode: {
      auto step_result = executor_.JoinANode(old_table,id_caches,one_step.join_node_);
      step_table = get<1>(step_result).values_;
      break;
    }
    case StepOperation::JoinType::JoinTable: {
      throw string("Optimizer::DepthSearchOneLayer not suitable for JoinTable");
    }
    case StepOperation::JoinType::GenerateCandidates : {
      throw string("Optimizer::DepthSearchOneLayer not suitable for GenerateCandidates");
    }
    case StepOperation::JoinType::EdgeCheck: {
      auto step_result = executor_.ANodeEdgesConstraintFilter(one_step.edge_filter_, old_table,id_caches);
      step_table = get<1>(step_result).values_;
      break;
    }
    case StepOperation::JoinType::JoinTwoNodes: {
      auto step_result = executor_.JoinTwoNode(one_step.join_two_node_, old_table,id_caches);
      step_table = get<1>(step_result).values_;
      break;
    }
  }


  if(step_table->empty())
    return make_tuple(false, layer_result);

  /* deep in bottom */
  if( (unsigned )layer_count + 1 == query_plan->join_order_->size()) {
    result_number_till_now += step_table->size();
    return make_tuple(false, layer_result);
  }

  auto all_result = make_shared<TableContent>();
  layer_result.values_ = all_result;

  /*  go deeper by filling a node */
  for(const auto& one_result:*step_table) {
    auto one_record_table = make_shared<TableContent>();
    one_record_table->push_back(one_result);
    auto next_layer_result  =this->DepthSearchOneLayer(query_plan, layer_count + 1,result_number_till_now,  limit_number,
                                                       one_record_table, id_caches,table_template);

    /* if not success */
    if(!get<0>(next_layer_result))
      continue;

    auto record_extended = get<1>(next_layer_result).values_;
    if(record_extended->empty())
      continue;

    /* record the result */
    for(const auto& next_layer_one_result: *record_extended) {
      all_result->push_back(next_layer_one_result);
    }
    if(limit_number!=-1) {
      /* examine the number*/
      if (result_number_till_now >= limit_number)
        break;
    }
  }
  return make_tuple(true,layer_result);
}


tuple<bool,shared_ptr<IntermediateResult>> Optimizer::DoQuery(SPARQLquery &sparql_query,QueryInfo query_info) {
  return MergeBasicQuery(sparql_query);
}

tuple<bool, shared_ptr<IntermediateResult>> Optimizer::DoQuery(std::shared_ptr<BGPQuery> bgp_query,QueryInfo query_info,
                                                               TopKStrategy top_k_method) {

#ifdef TOPK_DEBUG_INFO
  std::cout<<"Optimizer:: limit used:"<<query_info.limit_<<std::endl;
  std::cout<<"Optimizer::DoQuery limit num:"<<query_info.limit_num_<<std::endl;
#endif


  auto var_candidates_cache = make_shared<unordered_map<TYPE_ENTITY_LITERAL_ID,shared_ptr<IDList>>>();
  shared_ptr<QueryPlan> query_plan;
  auto strategy = this->ChooseStrategy(bgp_query,&query_info);
  if(strategy == BasicQueryStrategy::Normal)
  {

    long t1 =Util::get_cur_time();
    auto const_candidates = QueryPlan::OnlyConstFilter(bgp_query, this->kv_store_);
    for (auto &constant_generating_step: *const_candidates) {
      executor_.CacheConstantCandidates(constant_generating_step, var_candidates_cache);
    };
    long t2 = Util::get_cur_time();
    cout << "get var cache, used " << (t2 - t1) << "ms." << endl;

    long t3 = Util::get_cur_time();

    cout << "id_list.size = " << var_candidates_cache->size() << endl;
    //for (auto x : *var_candidates_cache){
    //  cout << "var[" << x.first << "] = " << bgp_query->get_vardescrip_by_id(x.first). << ", var_candidate list size = " << x.second->size() << endl;
    //}
#ifdef FEED_PLAN
    vector<int> node_order = {2,1,0};
      auto best_plan_tree = new PlanTree(node_order);
#else
    auto best_plan_tree = (new PlanGenerator(kv_store_, bgp_query.get(), statistics, var_candidates_cache))->get_random_plan();
    // get_plan(basic_query_pointer, this->kv_store_, var_candidates_cache);
#endif
    long t4 = Util::get_cur_time();
    cout << "plan get, used " << (t4 - t3) << "ms." << endl;
    best_plan_tree->print(bgp_query.get());
    cout << "plan print done" << endl;

    auto bfs_result = this->ExecutionBreathFirst(bgp_query,query_info,best_plan_tree->root_node,var_candidates_cache);

    long t5 = Util::get_cur_time();
    cout << "execution, used " << (t5 - t3) << "ms." << endl;

    auto bfs_table = get<1>(bfs_result);
    auto pos_var_mapping = bfs_table.pos_id_map;
    auto var_pos_mapping = bfs_table.id_pos_map;

    long t6 = Util::get_cur_time();
    CopyToResult(bgp_query->get_result_list_pointer(), bgp_query, bfs_table);
#ifdef OPTIMIZER_DEBUG_INFO
    cout<<"after copy bfs result size "<<bgp_query->get_result_list_pointer()->size()<<endl;
#endif
    long t7 = Util::get_cur_time();
    cout << "copy to result, used " << (t7 - t6) <<"ms." <<endl;
    cout << "total execution, used " << (t7 - t1) <<"ms."<<endl;
  }
  else if(strategy ==BasicQueryStrategy::Special){
    // if(bgp_query->get_triple_num() == 1 && bgp_query->get_total_var_num() == 3){
      // todo: get all triples in database
      //;

    printf("BasicQueryStrategy::Special not supported yet\n");
  }
  else if(strategy == BasicQueryStrategy::TopK)
  {
      long t_begin = Util::get_cur_time();
#ifdef TOPK_SUPPORT
    auto const_candidates = QueryPlan::OnlyConstFilter(bgp_query, this->kv_store_);
    for (auto &constant_generating_step: *const_candidates) {
      executor_.CacheConstantCandidates(constant_generating_step, var_candidates_cache);
    };

#ifdef TOPK_DEBUG_INFO
    std::cout<<"Top-k Constant Filtering Candidates Info"<<std::endl;
    for(const auto& pair:*var_candidates_cache)
      std::cout<<"var["<<pair.first<<"]  "<<pair.second->size()<<std::endl;
    std::cout<<"Top-k Constant Filtering Candidates Info End"<<std::endl;
#endif

    auto search_plan = make_shared<TopKSearchPlan>(bgp_query, this->kv_store_,
                                                   this->statistics, (*query_info.ordered_by_expressions_)[0],
                                                   var_candidates_cache);

#ifdef TOPK_DEBUG_INFO
    std::cout<<"Top-k Search Plan"<<std::endl;
    search_plan->DebugInfo(bgp_query,kv_store_);
#endif

    if(k_value_feed)
      query_info.limit_num_ = benchmark_k;
    tuple<bool, IntermediateResult> top_k_result;

    auto first_item = (*query_info.ordered_by_expressions_)[0];
    auto var_coefficients = TopKUtil::getVarCoefficients(first_item);
    if (first_item.descending)
      std::for_each(var_coefficients->begin(),
                    var_coefficients->end(),
                    [](decltype(*var_coefficients->begin()) pair_it) {
                      pair_it.second = -pair_it.second;
                    });

    search_plan->GetPlan(bgp_query, this->kv_store_,
                         this->statistics, (*query_info.ordered_by_expressions_)[0],
                         var_candidates_cache,var_coefficients);
    //search_plan->GetPlanCostModel(bgp_query, this->kv_store_,
    //                              this->statistics, (*query_info.ordered_by_expressions_)[0],
    //                              var_candidates_cache,query_info.limit_num_,var_coefficients);

    if(top_k_method==TopKStrategy::DP_B)
      top_k_result = this->ExecutionTopK_B(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::DP_P)
      top_k_result = this->ExecutionTopK_P(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::K_TPM)
      top_k_result = this->ExecutionTopK_kTPM(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::SAE)
      top_k_result = this->ExecutionTopK_SAE(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::DFS_OLD) // no compress vector of DFS
      top_k_result = this->ExecutionTopK_DfsOld(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::PTAB)
      top_k_result = this->ExecutionTopK_Dfs(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::TAKE_ALL)
      top_k_result = this->ExecutionTopK_TakeAll(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::TAKE2)
      top_k_result = this->ExecutionTopK_Take2(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method==TopKStrategy::EAGER)
       top_k_result = this->ExecutionTopK_Eager(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    if(top_k_method == TopKStrategy::RankAfterMatching)
      top_k_result = this->ExecutionTopK_RAM(bgp_query, search_plan, query_info, var_candidates_cache,var_coefficients);
    auto result_table = get<1>(top_k_result);

    long t_end = Util::get_cur_time();
    cout << "TOPK execution, used " << t_end - t_begin << "ms." << endl;
    CopyToResult(bgp_query->get_result_list_pointer(), bgp_query, result_table);
#endif
  }
  else if(strategy == BasicQueryStrategy::limitK)
  {

  }
  return tuple<bool, shared_ptr<IntermediateResult>>();

}

tuple<bool, shared_ptr<IntermediateResult>> Optimizer::MergeBasicQuery(SPARQLquery &sparql_query) {
  return tuple<bool, shared_ptr<IntermediateResult>>();
}

/**
 * 调用栈：   Strategy::handle
 *          Strategy::handler0
 *          Join::join_basic
 * @param target
 * @param basic_query
 * @param result, which has already gotten CoreVar, Pre var(selected), but not get satellites
 * @return
 */
bool Optimizer::CopyToResult(vector<unsigned int *> *target,
                             BasicQuery* basic_query,
                             const shared_ptr<IntermediateResult>& result) {

  assert(target->empty());

  int select_var_num = basic_query->getSelectVarNum();
  int core_var_num = basic_query->getRetrievedVarNum();
  int selected_pre_var_num = basic_query->getSelectedPreVarNum();

  // A little different with that in Join::CopyToResult
  // may basic query don't allocate an id for not selected var so selected_pre_var_num = pre_var_num?
  cout << "position to var des size: " << result->pos_id_map->size() << endl;
  cout << "totalvar_num: " << basic_query->getTotalVarNum()<<endl;
  cout << "varnum: " << basic_query->getVarNum()<<endl;
  cout << "selected var num: " << select_var_num<<endl;
  cout << "core var num: " << core_var_num<<endl;
  cout << "selected pre var num: " << selected_pre_var_num<<endl;

#ifdef OPTIMIZER_DEBUG_INFO
  for(auto &pos_id:*result->pos_id_map)
  {
    cout<<"pos["<<pos_id.first<<"] "<<pos_id.second<<endl;
  }
#endif
//  Linglin Yang fix it to total_var_num,
//  maybe change it to core_var_num + selected_pre_var_num in the future
  if (result->pos_id_map->size() != basic_query->getVarNum())
  {
    cout << "terrible error in Optimizer::CopyToResult!" << endl;
    return false;
  }

  shared_ptr<vector<Satellite>> satellites = make_shared<vector<Satellite>>();

  auto record_len = select_var_num + selected_pre_var_num;
  auto record = new unsigned[record_len];

  auto position_id_map_ptr = result->pos_id_map;
  auto id_position_map_ptr = result->id_pos_map;
  cout<<"fir_var_position"<<basic_query->getSelectedVarPosition((*position_id_map_ptr)[0])<<endl;
  int var_num =  basic_query->getVarNum();
  for (const auto&  record_ptr : *(result->values_))
  {
    int column_index = 0;
    for (; column_index < core_var_num; ++column_index)
    {
      //This is because selected var id is always smaller
      if ( (*position_id_map_ptr)[column_index] < select_var_num)
      {
        int vpos = basic_query->getSelectedVarPosition((*position_id_map_ptr)[column_index]);
        record[vpos] = (*record_ptr)[column_index];
      }
    }

    // below are for selected pre vars
    while (column_index < result->pos_id_map->size() )
    {
      //only add selected ones
      // 原句是 int pre_var_id = this->pos2id[i] - this->var_num;
      // 可能有bug
      int pre_var_id = (*position_id_map_ptr)[column_index] - var_num;
      int pre_var_position = basic_query->getSelectedPreVarPosition(pre_var_id);
      if(pre_var_position >= 0)
      {
        record[pre_var_position] = (*record_ptr)[column_index];
      }
      ++column_index;
    }

    bool valid = true;
    //generate satellites when constructing records
    //NOTICE: satellites in join must be selected
    //core vertex maybe not in select
    for (column_index = 0; column_index < core_var_num; ++column_index)
    {
      int id = (*position_id_map_ptr)[column_index];
      unsigned ele = (*record_ptr)[column_index];
      int degree = basic_query->getVarDegree(id);
      for (int j = 0; j < degree; ++j)
      {
        int id2 = basic_query->getEdgeNeighborID(id, j);
        if (basic_query->isSatelliteInJoin(id2) == false)
          continue;

        if(id_position_map_ptr->find(id2)!=id_position_map_ptr->end())
          continue;

        cout << "satellite node: " <<basic_query->getVarName(id2)<<endl;
        unsigned* idlist = nullptr;
        unsigned idlist_len = 0;
        int triple_id = basic_query->getEdgeID(id, j);
        Triple triple = basic_query->getTriple(triple_id);

        TYPE_PREDICATE_ID preid = basic_query->getEdgePreID(id, j);
        if (preid == -2)  //?p
        {
          string predicate = triple.predicate;
          int pre_var_id = basic_query->getPreVarID(predicate);
          //if(this->basic_query->isPreVarSelected(pre_var_id))
          //{
          preid = (*record_ptr)[ (*id_position_map_ptr)[pre_var_id+ var_num]];
          //}
        }
        else if (preid == -1)  //INVALID_PREDICATE_ID
        {
          //ERROR
        }

        char edge_type = basic_query->getEdgeType(id, j);
        if (edge_type == Util::EDGE_OUT)
        {
          this->kv_store_->getobjIDlistBysubIDpreID(ele, preid, idlist, idlist_len, true, this->txn_);
        }
        else
        {
          this->kv_store_->getsubIDlistByobjIDpreID(ele, preid, idlist, idlist_len, true, this->txn_);
        }

        if(idlist_len == 0)
        {
          valid = false;
          break;
        }
        satellites->push_back(Satellite(id2, idlist, idlist_len));

      }
      if(!valid)
      {
        break;
      }
    }

    int size = satellites->size();
    if(valid)
    {
      // Join::cartesian
      Cartesian(0,size,record_len,record,satellites,target,basic_query);
    }

    for (int k = 0; k < size; ++k)
    {
      delete[] (*satellites)[k].idlist;
      //this->satellites[k].idlist = NULL;
    }
    //WARN:use this to avoid influence on the next loop
    satellites->clear();
  }

  delete[] record;


  return true;
}


bool Optimizer::CopyToResult(vector<unsigned int *> *target,
                             shared_ptr<BGPQuery> bgp_query,
                             IntermediateResult result) {
  return true;
  assert(target->empty());

#ifdef OPTIMIZER_DEBUG_INFO
  cout << "position to var des size: " << result.pos_id_map->size() << endl;
  cout << "total var_num: " << bgp_query->get_total_var_num()<<endl;
  cout << "selected var num: " << bgp_query->selected_so_var_num<<endl;
  cout<<"Before Copy, result size:"<<result.values_->size()<<endl;
  for(auto &pos_id:*result.pos_id_map)
    cout<<"pos["<<pos_id.first<<"] "<<pos_id.second<<"\t";
  cout<<endl;
#endif
  auto record_len = bgp_query->selected_so_var_num + bgp_query->selected_pre_var_num;
  auto record = new unsigned[record_len];

  auto position_id_map_ptr = result.pos_id_map;
  auto id_position_map_ptr = result.id_pos_map;

  // position_map[i] means in the new table, the i-th column
  // is the  position_map[i]-th column from old table
  auto position_map = new TYPE_ENTITY_LITERAL_ID[record_len];

  for(auto pos_id_pair:bgp_query->position_id_map)
  {
    auto var_id = pos_id_pair.second;
    auto old_position = (*id_position_map_ptr)[var_id];
    auto new_position =  pos_id_pair.first;
    position_map[new_position] = old_position;
  }

  for (const auto&  record_ptr : *(result.values_))
  {
    auto new_record = new unsigned[record_len];
    for (int column_index = 0; column_index < record_len; ++column_index)
    {
      auto old_position = position_map[column_index];
      record[column_index] = (*record_ptr)[old_position];
    }
    target->push_back(new_record);
  }

  delete[] position_map;
  return false;
}


void
Optimizer::Cartesian(int pos, int end,int record_len,unsigned* record,
                     const shared_ptr<vector<Satellite>>& satellites,
                     vector<unsigned*>* result_list,
                     BasicQuery *basic_query)
{
  if (pos == end)
  {
    auto new_record = new unsigned[record_len];
    memcpy(new_record, record, sizeof(unsigned) * record_len);
    result_list->push_back(new_record);
    return;
  }

  unsigned size = (*satellites)[pos].idlist_len;
  int id = (*satellites)[pos].id;
  int vpos = basic_query->getSelectedVarPosition(id);
  unsigned* list = (*satellites)[pos].idlist;
  for (unsigned i = 0; i < size; ++i)
  {
    record[vpos] = list[i];
    Cartesian(pos + 1, end,record_len,record,satellites,result_list,basic_query);
  }
}

/**
 * It is not implemented
 */
tuple<bool,PositionValueSharedPtr, TableContentShardPtr> Optimizer::ExecutionBreathFirst(BasicQuery* basic_query,
                                                                                         const QueryInfo& query_info,
                                                                                         Tree_node* plan_tree_node,
                                                                                         IDCachesSharePtr id_caches)
{
  return make_tuple(false, nullptr,nullptr);
}

/**
 *
 * @param basic_query
 * @param query_info
 * @param plan_tree_node
 * @param id_caches
 * @return PositionValueSharedPtr position2var
 */
tuple<bool,IntermediateResult> Optimizer::ExecutionBreathFirst(shared_ptr<BGPQuery> bgp_query,
                                                               QueryInfo query_info,
                                                               Tree_node* plan_tree_node,
                                                               IDCachesSharePtr id_caches)
{
  // leaf node
  auto step_operation = plan_tree_node->node;
  auto operation_type = step_operation->join_type_;

  // leaf node : create a table
  if (plan_tree_node->left_node == nullptr && plan_tree_node->right_node == nullptr)
  {
    auto position_id_map = make_shared<PositionValue>();
    IntermediateResult leaf_table;
    if(operation_type== StepOperation::JoinType::JoinNode){
      auto r = this->PrepareInitial(bgp_query,step_operation->join_node_);
      bool is_entity= get<0>(r);
      bool is_literal = get<1>(r);
      bool is_predicate = get<2>(r);
      auto initial_result = executor_.InitialTableOneNode(step_operation->join_node_,is_entity,is_literal,is_predicate,id_caches);
      leaf_table = get<1>(initial_result);
#ifdef OPTIMIZER_DEBUG_INFO
      cout<<"JoinNode result size:"<<leaf_table.values_->size()<<endl;
#endif
    }
    else if(operation_type== StepOperation::JoinType::GenerateCandidates)
    {
      executor_.UpdateIDCache(step_operation,id_caches);
    }
    else if(operation_type==StepOperation::JoinType::JoinTable)
      throw string("StepOperation::JoinType::JoinTable cannot happened in leaf node");
    else if(operation_type==StepOperation::JoinType::EdgeCheck)
      throw string("StepOperation::JoinType::EdgeCheck cannot happened in leaf node");
    else if(operation_type==StepOperation::JoinType::JoinTwoNodes){
      auto initial_result = executor_.InitialTableTwoNode(step_operation->join_two_node_,id_caches);
      leaf_table = get<1>(initial_result);
    }
    else
      throw string("unexpected JoinType");

    return make_tuple(true,leaf_table);
  }


  //inner node
  tuple<bool,IntermediateResult> left_r;
  tuple<bool,IntermediateResult> right_r;

  if(plan_tree_node->left_node != nullptr)
    left_r = this->ExecutionBreathFirst(bgp_query,query_info,plan_tree_node->left_node,id_caches);

  if(plan_tree_node->right_node != nullptr)
    right_r = this->ExecutionBreathFirst(bgp_query,query_info,plan_tree_node->right_node,id_caches);

  if(operation_type == StepOperation::JoinType::GenerateCandidates){
    executor_.UpdateIDCache(step_operation,id_caches);
    return left_r;
  }

  auto left_table = get<1>(left_r);
  // reading the left child result
  auto left_records = left_table.values_;

  // JoinNode/EdgeCheck/JoinTwoNodes has only ONE left child
  if(operation_type== StepOperation::JoinType::JoinNode) {
    // create a new table
    if(left_table.GetColumns()==0){
      auto position_id_map = make_shared<PositionValue>();
      auto r = this->PrepareInitial(bgp_query,step_operation->join_node_);
      bool is_entity= get<0>(r);
      bool is_literal = get<1>(r);
      bool is_predicate = get<2>(r);
      auto initial_result = executor_.InitialTableOneNode(step_operation->join_node_,is_entity,is_literal,is_predicate,id_caches);
      auto leaf_table = get<1>(initial_result);
      return make_tuple(true,leaf_table);
    }

    auto one_step_join = plan_tree_node->node;
    auto node_to_join = one_step_join->join_node_->node_to_join_;
    cout<<"join node ["<<bgp_query->get_var_name_by_id(node_to_join)<<"]"<<",  ";
    auto join_node = one_step_join->join_node_;
    auto edges = *join_node->edges_;
    auto edge_c = *join_node->edges_constant_info_;
#ifdef TABLE_OPERATOR_DEBUG_INFO
    for(int i=0;i<edges.size();i++)
      {
        std::cout<<"edge["<<i<<"] "<<edges[i].toString()<<std::endl;
        std::cout<<"constant["<<i<<"] "<<edge_c[i].toString()<<std::endl;
      }
#endif
    long t1 = Util::get_cur_time();
    auto step_result = executor_.JoinANode(left_table,id_caches,one_step_join->join_node_);

#ifdef OPTIMIZER_DEBUG_INFO
    cout<<"JoinNode result size:"<<get<1>(step_result).values_->size()<<endl;
#endif
    long t2 = Util::get_cur_time();
    cout<< ",  used " << (t2 - t1) << "ms." <<endl;
    return step_result;
  }
  else if(operation_type==StepOperation::JoinType::EdgeCheck){

    auto edge_filter = step_operation->edge_filter_;
    long t1 = Util::get_cur_time();
    auto step_result =executor_.ANodeEdgesConstraintFilter(edge_filter, left_table, id_caches);

    cout<<"result size "<<get<1>(step_result).values_->size();
    long t2 = Util::get_cur_time();
    cout<< ",  used " << (t2 - t1) << "ms." <<endl;
    return step_result;
  }
  else if(operation_type==StepOperation::JoinType::JoinTwoNodes){
    if(left_table.GetColumns()==0){
      auto position_id_map = make_shared<PositionValue>();
      auto initial_result = executor_.InitialTableTwoNode(step_operation->join_two_node_,id_caches);
#ifdef OPTIMIZER_DEBUG_INFO
      cout<<"JoinTwoNodes result size:"<<get<1>(initial_result).values_->size()<<endl;
#endif
      return initial_result;
    }
    auto one_step_join = plan_tree_node->node;
    auto join_two_plan = one_step_join->join_two_node_;

    auto node1 = join_two_plan->node_to_join_1_;
    auto node2 = join_two_plan->node_to_join_2_;
    cout<<"join node ["<<bgp_query->get_var_name_by_id(node1)<<"]"<<",  ";
    cout<<"join node ["<<bgp_query->get_var_name_by_id(node2)<<"]"<<",  ";

    long t1 = Util::get_cur_time();
    auto step_result = executor_.JoinTwoNode(join_two_plan, left_table, id_caches);

    long t2 = Util::get_cur_time();
    cout<< ",  used " << (t2 - t1) << "ms." <<endl;

#ifdef OPTIMIZER_DEBUG_INFO
    cout<<"JoinTwoNodes result size:"<<get<1>(step_result).values_->size()<<endl;
#endif
    return step_result;
  }

  auto right_table = get<1>(right_r);
  if(operation_type==StepOperation::JoinType::JoinTable){
    return executor_.JoinTable(step_operation->join_table_, left_table, right_table);
  }
  else
    throw string("unexpected JoinType");
  return make_tuple(false,left_table);
}

/**
 * calculator which role a variable node maybe
 * @param bgp_query
 * @param join_a_node_plan
 * @return is_entity,is_literal,is_predicate
 */
tuple<bool,bool,bool>
Optimizer::PrepareInitial(shared_ptr<BGPQuery> bgp_query,
                          shared_ptr<FeedOneNode> join_a_node_plan) const {
  bool is_entity= false;
  bool is_predicate= false;
  bool is_literal= false;
  auto target_var_id = join_a_node_plan->node_to_join_;
  cout << "leaf node [" << bgp_query->get_var_name_by_id(target_var_id) << "],  ";
  auto var_descriptor = bgp_query->get_vardescrip_by_id(target_var_id);
  if (var_descriptor->var_type_ == VarDescriptor::VarType::Predicate) {
    is_predicate = true;
  }
  else {
    is_entity = true;
    auto var_name = bgp_query->get_var_name_by_id(target_var_id);
    auto edge_ids = var_descriptor->so_edge_index_;
    for (auto edge_id : edge_ids) {
      auto triple = bgp_query->get_triple_by_index(edge_id);
      if (var_name == triple.getObject()) {
        is_literal = true;
        break;
      }
    }
  }
  return make_tuple(is_entity,is_literal,is_predicate);
}

#ifdef TOPK_SUPPORT
tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_B(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                          const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                          std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients) {

  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store = this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout << " Optimizer::ExecutionTopK_B  env->id_caches " << env->id_caches->size() << endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();

  test_scores.clear();
  result_ids.clear();
#ifdef ITERATOR_COUNT
  FQ_NUM=0;
  FR_NUM=0;
  OW_NUM=0;
#endif
  auto t0 = Util::get_cur_time();
  auto root_fr = DPBUtil::BuildIteratorTree(tree_search_plan, env);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for (auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  for (int i = 1; i <= query_info.limit_num_ + deleted_num; i++) {
    root_fr->TryGetNext(k + deleted_num);
    if (root_fr->pool_.size() != (unsigned int) i)
      break;
#ifdef TOPK_DEBUG_INFO
    else {
      cout << "get top-" << i << " "<<root_fr->pool_[i-1].cost<<endl;
    }
#endif

    // can't get any more
    if ((unsigned) i > root_fr->pool_.size())
      break;

    auto var_num = pos_var_mapping->size();

    auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    record->reserve(var_num);
    root_fr->GetResult(i - 1, record);

    bool success = true;

//    if(tree_search_plan->HasCycle())
//    {
//      auto &non_tree_edge = tree_search_plan->GetNonTreeEdges();
//      auto temp_content_ptr = make_shared<TableContent>();
//      temp_content_ptr->push_back(record);
//      temp_table.values_ = temp_content_ptr;
//      auto filter_result = this->executor_.ANodeEdgesConstraintFilter(non_tree_edge.edge_filter_,temp_table,id_caches);
//      auto filter_result_size = get<1>(filter_result).values_->size();
//      success = filter_result_size != 0;
//    }

    if (success) {
      result_list->push_back(record);
      if(any_k && result_list->size()%any_k_step==0)
        any_k_time.push_back(Util::get_cur_time()-t0);
    }
    else
      deleted_num += 1;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1 - t0;

#ifdef TOPK_DEBUG_RESULT_INFO
  auto var_num = bgp_query->get_total_var_num();
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(unsigned i =0;i<result_list->size(); i++) {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<root_fr->pool_[i].cost;
    test_scores.push_back(root_fr->pool_[i].cost);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    // cout<<endl;
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<(*rec)[j];
      else
        cout<<" "<<(*rec)[j];
    cout<<endl;
    it++;
  }
#endif
#ifdef ITERATOR_COUNT
  if (iterator_count_enable) {
    fq_num_vec.push_back(FQ_NUM);
    fr_num_vec.push_back(FR_NUM);
    ow_num_vec.push_back(OW_NUM);
  }
#endif
  delete root_fr;
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_P(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                          const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                          std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_P  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  // create an empty global queue (min-heap), GQ;
  GlobalQueue global_queue;
  // corresponding to line 2-3 in Algo.5, that is,
  // for each v.OWJ (v is a data node in GR) do
  //     e = v.OWJ.get next() ; if e != END then GQ.insert(e);
  auto t0 = Util::get_cur_time();
  auto root_fr = DPPUtil::BuildIteratorTree(tree_search_plan, env,&global_queue);
#ifdef TOPK_DEBUG_INFO
  vector<unsigned int> id_occur_times(env->bgp_query->get_total_var_num(),0);
  cout<<"subtree_has_coefficient ";
  for(int i = 0;i<env->bgp_query->get_total_var_num();i++)
    cout<<" "<<i<<"["<<env->subtree_has_coefficient_[i]<<"]\t";
  cout<<endl;
#endif
  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  vector<double> scores_list_;
  // corresponding to line 5 in Algo.5
  while(!global_queue.Empty() && result_list->size()<k)
  {
    auto e = global_queue.ExtractMin();
    // corresponding to line 7 in Algo.5
    if(e->type_ == QueueElement::Type::FR || e->type_ == QueueElement::Type::OW)
    {
      // corresponding to line 7 in Algo.5 , that is,
      //    v.iPoolJ.size++; v.iPoolJ[v.iPoolJ.size] = e;
      //    trigger_seq(v, J);
      auto parent_fq = dynamic_cast<DPPFQIterator*>( e->fq_iterator_);
      auto self_fr_or_ow = e->self_dpp_list_;
      if(e->type_ == QueueElement::Type::FR)
        self_fr_or_ow->pool_.push_back(*(e->i_element_));
#ifdef TOPK_DEBUG_INFO
      if(e->type_ == QueueElement::Type::OW)
        id_occur_times[self_fr_or_ow->var_id]++;
#endif
      DPPUtil::TriggerSeq(env->k,parent_fq,self_fr_or_ow,e->i_th_type,&global_queue);
    }
    else  // FQ
    {
      // corresponding to line 10 in Algo.5
      //v.ePool.size++; v.ePool[v.ePool.size] = e;
      auto fq_itself = dynamic_cast<DPPFQIterator*>( e->fq_iterator_);
      fq_itself->InsertEPool(*e->fq_element_);
#ifdef TOPK_DEBUG_INFO
      id_occur_times[fq_itself->var_id]++;
#endif
      // corresponding to line 12 in Algo.5
      if(fq_itself->var_id == root_var)
      {
        // corresponding to line 13-14 in Algo.5
        //  i++; derive Ai by backtracking e;
        //  if i < k then v.FQ.get next() ; else break
        auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
        record->reserve(var_num);
        auto already_result_num = fq_itself->EPoolSize();
        fq_itself->GetResult(already_result_num-1,record);
        result_list->push_back(record);
        scores_list_.push_back(e->Score());
        if(result_list->size()<k)
          fq_itself->TryGetNext(k,&global_queue);
      }
      else     // corresponding to line 15 in Algo.5
        DPPUtil::RorDExpansion(env->k,fq_itself,&global_queue,tree_search_plan.get());
    }
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];

    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }
#endif
  delete root_fr;
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_kTPM(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                          const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                             std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_kTPM  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = kTPMUtil::BuildIteratorTree(tree_search_plan, env);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  vector<double> scores_list_;
  auto root_space = kTPMUtil::WrapStartSpace(root_fr);
  SpaceHeap space_heap;
  space_heap.push(std::move(root_space));
  while(!space_heap.empty())
  {
    auto best_space = space_heap.popMin();
    scores_list_.push_back(best_space->best_score_);
    auto record = kTPMUtil::Divide(k,best_space,&space_heap);
    result_list->push_back(record);
    Util::get_cur_time() - t0;
    if(any_k && result_list->size()%any_k_step==0)
      any_k_time.push_back(Util::get_cur_time()-t0);
    if(result_list->size()>=(unsigned int)k)
      break;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];
    test_scores.push_back(scores_list_[i]);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }
#endif
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_SAE(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                             const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                             std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_SAE  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = SAEUtil::BuildIteratorTree(tree_search_plan, env);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  vector<double> scores_list_;
  auto root_space = SAEUtil::WrapStartSpace(root_fr);
  SpaceHeapSAE space_heap;
  space_heap.push(std::move(root_space));
  while(!space_heap.empty())
  {
    auto best_space = space_heap.popMin();
    scores_list_.push_back(best_space->best_score_);
    auto record = SAEUtil::Divide(k,best_space,&space_heap);
    bool success  = true;
    if(tree_search_plan->HasCycle())
    {
      auto &non_tree_edge = tree_search_plan->GetNonTreeEdges();
      auto temp_content_ptr = make_shared<TableContent>();
      temp_content_ptr->push_back(record);
      temp_table.values_ = temp_content_ptr;
      auto filter_result = this->executor_.ANodeEdgesConstraintFilter(non_tree_edge.edge_filter_,temp_table,id_caches);
      auto filter_result_size = get<1>(filter_result).values_->size();
      success = filter_result_size != 0;
    }
    if(success){
      result_list->push_back(record);
      Util::get_cur_time() - t0;
      if (any_k && result_list->size() % any_k_step == 0)
        any_k_time.push_back(Util::get_cur_time() - t0);
      if (result_list->size() >= (unsigned int) k)
        break;
    }
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];
    test_scores.push_back(scores_list_[i]);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }
#endif
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_DfsOld(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                            const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                               std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients)
{
  // CompressedVector::InitialCombinatorial(2040,8);
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_DfsOld  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  auto dfs_help_info = new DfsHelpInfo;
  NodeCaches node_caches;
  node_caches.var_valid_ids = std::unique_ptr<std::unordered_map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DfsList>>[]>(
      new std::unordered_map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DfsList>>[bgp_query->get_total_var_num()]);
  node_caches.var_invalid_ids = std::unique_ptr<std::unordered_set<TYPE_ENTITY_LITERAL_ID>[]>(
      new std::unordered_set<TYPE_ENTITY_LITERAL_ID>[bgp_query->get_total_var_num()]);
  dfs_help_info->node_caches = &node_caches;
  dfs_help_info->hop_index = hop_index_;
  std::unique_ptr<TopKTreeNode*> node_plan_list =  std::unique_ptr<TopKTreeNode*>(
      new TopKTreeNode*[bgp_query->get_total_var_num()]
  );
  TopKUtil::FillArrayOfTreeNode(tree_search_plan->tree_root_,node_plan_list.get());
  dfs_help_info->node_plan = node_plan_list.get();

  dfs_help_info->env = env;
  dfs_help_info->kv_store = env->kv_store;
  dfs_help_info->basic_query = env->basic_query;
  dfs_help_info->bgp_query = env->bgp_query;
  dfs_help_info->id_caches = env->id_caches;
  dfs_help_info->k = env->k;
  dfs_help_info->non_tree_edges_lists_ = env->non_tree_edges_lists_;
  dfs_help_info->coefficients = env->coefficients;
  dfs_help_info->txn = env->txn;
  dfs_help_info->ss = env->ss;
  dfs_help_info->ow_coefficient_counts_ = env->ow_coefficient_counts_;
  dfs_help_info->max_leaf_distance_ = tree_search_plan->GetMaxLeafDistance();
  dfs_help_info->DfsNew = false;
  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = DfsUtilDynamic::BuildIteratorTree(tree_search_plan, dfs_help_info);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif
  for(int i =1; i<=query_info.limit_num_ + deleted_num; i++)
  {
    root_fr->TryGetNext(dfs_help_info);
    if(root_fr->pool_.size()!=(unsigned int)i)
      break;
#ifdef TOPK_DEBUG_INFO
    else {
      cout << "get top-" << i << " "<<root_fr->pool_[i-1].cost<<endl;
    }
#endif

    // can't get any more
    if((unsigned)i>root_fr->pool_.size())
      break;

    auto var_num = pos_var_mapping->size();

    auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    record->reserve(var_num);
    root_fr->GetResult(i-1,record);

    bool success = true;

    if(success){
      result_list->push_back(record);
      if(any_k && result_list->size()%any_k_step==0)
        any_k_time.push_back(Util::get_cur_time()-t0);
    }
    else
      deleted_num += 1;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<root_fr->pool_[i].cost;
    test_scores.push_back(root_fr->pool_[i].cost);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    // cout<<endl;
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<(*rec)[j];
      else
        cout<<" "<<(*rec)[j];
    cout<<endl;
    it++;
  }
#endif
  delete env;
  delete dfs_help_info;
  return std::make_tuple(true,result_table);
}

tuple<bool, IntermediateResult>
Optimizer::ExecutionTopK_Dfs(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                             const QueryInfo &query_info, IDCachesSharePtr id_caches,
                             std::shared_ptr<std::unordered_map<std::string, double>> var_coefficients) {
  CompressedVector::InitialCombinatorial(1000,10);
#ifdef ITERATOR_COUNT
  FQ_NUM=0;
  FR_NUM=0;
  OW_NUM=0;
#endif
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout << " Optimizer::ExecutionTopK_Dfs  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  auto dfs_help_info = new DfsHelpInfo;
  NodeCaches node_caches;
  node_caches.var_valid_ids = std::unique_ptr<std::unordered_map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DfsList>>[]>(
      new std::unordered_map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DfsList>>[bgp_query->get_total_var_num()]);
  node_caches.var_invalid_ids = std::unique_ptr<std::unordered_set<TYPE_ENTITY_LITERAL_ID>[]>(
      new std::unordered_set<TYPE_ENTITY_LITERAL_ID>[bgp_query->get_total_var_num()]);
  dfs_help_info->node_caches = &node_caches;
  dfs_help_info->hop_index = hop_index_;
  std::unique_ptr<TopKTreeNode*> node_plan_list =  std::unique_ptr<TopKTreeNode*>(
      new TopKTreeNode*[bgp_query->get_total_var_num()]
      );
  TopKUtil::FillArrayOfTreeNode(tree_search_plan->tree_root_,node_plan_list.get());
  dfs_help_info->node_plan = node_plan_list.get();

  dfs_help_info->env = env;
  dfs_help_info->kv_store = env->kv_store;
  dfs_help_info->basic_query = env->basic_query;
  dfs_help_info->bgp_query = env->bgp_query;
  dfs_help_info->id_caches = env->id_caches;
  dfs_help_info->k = env->k;
  dfs_help_info->non_tree_edges_lists_ = env->non_tree_edges_lists_;
  dfs_help_info->coefficients = env->coefficients;
  dfs_help_info->txn = env->txn;
  dfs_help_info->ss = env->ss;
  dfs_help_info->ow_coefficient_counts_ = env->ow_coefficient_counts_;
  dfs_help_info->max_leaf_distance_ = tree_search_plan->GetMaxLeafDistance();
  dfs_help_info->DfsNew = true;
  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = DfsUtilCompressedVector::BuildIteratorTree(tree_search_plan, dfs_help_info);

  int deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();
#ifdef ITERATOR_COUNT
  if (iterator_count_enable) {
    fq_num_vec.push_back(FQ_NUM);
    fr_num_vec.push_back(FR_NUM);
    ow_num_vec.push_back(OW_NUM);
  }
#endif
#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif
  for (int i = 1; i <= query_info.limit_num_ + deleted_num; ++i) {
    root_fr->TryGetNext(dfs_help_info);
    if(root_fr->pool_.size()!=(unsigned int)i)
      break;
#ifdef TOPK_DEBUG_INFO
    else {
      cout << "get top-" << i << " "<<root_fr->pool_[i-1].cost<<endl;
    }
#endif

    // can't get anymore
    if((unsigned)i > root_fr->pool_.size())
      break;

    auto record = make_shared<vector<TYPE_ENTITY_LITERAL_ID>>();
    record->reserve(var_num);
    root_fr->GetResult(i-1,record);

    bool success = true;

    if (success) {
      result_list->push_back(record);
      if(any_k && result_list->size()%any_k_step==0)
        any_k_time.push_back(Util::get_cur_time()-t0);
#ifdef ITERATOR_COUNT
      if (iterator_count_enable) {
        fq_num_vec.push_back(FQ_NUM);
        fr_num_vec.push_back(FR_NUM);
        ow_num_vec.push_back(OW_NUM);
      }
#endif
    } else deleted_num += 1;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for (unsigned i = 0; i < result_list->size(); ++i) {
    auto rec = *it;
    cout << " record[" << i << "]" << " score:" << root_fr->pool_[i].cost;
    test_scores.push_back(root_fr->pool_[i].cost);
    result_ids.push_back(*rec);
    for(unsigned j = 0; j < var_num; j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    // cout<<endl;
    for(unsigned j = 0; j < var_num; ++j)
      if (is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<(*rec)[j];
      else
        cout<<" "<<(*rec)[j];
    cout<<endl;
    it++;
  }
#endif
  delete env;
  delete dfs_help_info;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_TakeAll(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                                const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                                std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;

  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_TakeAll  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = TakeAllUtil::BuildIteratorTree(tree_search_plan, env);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  vector<double> scores_list_;
  auto root_space = TakeAllUtil::WrapStartSpace(root_fr);
  TakeAllSpaceHeap space_heap;
  space_heap.push(std::move(root_space));
  while(!space_heap.empty())
  {
    auto best_space = space_heap.popMin();
    scores_list_.push_back(best_space->best_score_);
    auto record = TakeAllUtil::Divide(k,best_space,&space_heap);
    result_list->push_back(record);
    if(any_k && result_list->size()%any_k_step==0)
      any_k_time.push_back(Util::get_cur_time()-t0);
    if(result_list->size()>=(unsigned int)k)
      break;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];
    test_scores.push_back(scores_list_[i]);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }
#endif
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_Take2(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                              const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                              std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_Take2  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = Take2Util::BuildIteratorTree(tree_search_plan, env);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  vector<double> scores_list_;
  auto root_space = Take2Util::WrapStartSpace(root_fr);
  Take2SpaceHeap space_heap;
  space_heap.push(std::move(root_space));
  while(!space_heap.empty())
  {
    auto best_space = space_heap.popMin();
    scores_list_.push_back(best_space->best_score_);
    auto record = Take2Util::Divide(k,best_space,&space_heap);
    result_list->push_back(record);
    if(any_k && result_list->size()%any_k_step==0)
      any_k_time.push_back(Util::get_cur_time()-t0);
    if(result_list->size()>=(unsigned int)k)
      break;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];
    test_scores.push_back(scores_list_[i]);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }
#endif
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_Eager(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                              const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                              std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients){
  auto root_var = tree_search_plan->tree_root_->var_id;
  auto pos_var_mapping = TopKUtil::CalculatePosVarMapping(tree_search_plan);
  auto k = query_info.limit_num_;
  // Build Iterator tree
  auto env = new TopKUtil::Env();
  env->kv_store= this->kv_store_;
  env->bgp_query = bgp_query;
  env->id_caches = id_caches;
  cout<<" Optimizer::ExecutionTopK_Eager  env->id_caches "<<  env->id_caches->size()<<endl;
  env->k = query_info.limit_num_;
  env->coefficients = var_coefficients;
  env->txn = this->txn_;
  env->ss = make_shared<stringstream>();
  env->ow_coefficient_counts_.resize(bgp_query->get_total_var_num());
  env->subtree_has_coefficient_.resize(bgp_query->get_total_var_num(),false);

  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  auto root_fr = EagerUtil::BuildIteratorTree(tree_search_plan, env);

  decltype(query_info.limit_num_) deleted_num = 0;
  auto result_list = make_shared<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>();

  auto result_table = IntermediateResult();
  result_table.values_ = result_list;
  result_table.pos_id_map = pos_var_mapping;
  for(auto pos_id_pair:*pos_var_mapping)
    (*result_table.id_pos_map)[pos_id_pair.second] = pos_id_pair.first;

  auto temp_table = IntermediateResult();
  temp_table.pos_id_map = result_table.pos_id_map;
  temp_table.id_pos_map = result_table.id_pos_map;
  auto var_num = pos_var_mapping->size();

#ifdef TOPK_DEBUG_INFO
  for(const auto& pos_id_pair:*pos_var_mapping)
  {
    cout<<"pos["<<pos_id_pair.first<<"]"<<"|"<<pos_id_pair.second<<"| "<<bgp_query->get_var_name_by_id(pos_id_pair.second)<<"\t";
  }
  cout<<endl;
#endif

  vector<double> scores_list_;
  auto root_space = EagerUtil::WrapStartSpace(root_fr);
  EagerSpaceHeap space_heap;
  space_heap.push(std::move(root_space));
  while(!space_heap.empty())
  {
    auto best_space = space_heap.popMin();
    scores_list_.push_back(best_space->best_score_);
    auto record = EagerUtil::Divide(k,best_space,&space_heap);
    result_list->push_back(record);
    if(result_list->size()>=(unsigned int)k)
      break;
  }
  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;
#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];
    test_scores.push_back(scores_list_[i]);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*pos_var_mapping)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }
#endif
  delete env;
  return std::make_tuple(true,result_table);
}

tuple<bool,IntermediateResult> Optimizer::ExecutionTopK_RAM(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                              const QueryInfo& query_info, IDCachesSharePtr var_candidates_cache,
                                                              std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients)
{
  test_scores.clear();
  result_ids.clear();

  auto t0 = Util::get_cur_time();
  IntermediateResult r;

  {
    auto const_candidates = QueryPlan::OnlyConstFilter(bgp_query, this->kv_store_);
    for (auto &constant_generating_step: *const_candidates) {
      executor_.CacheConstantCandidates(constant_generating_step, var_candidates_cache);
    };


#ifdef FEED_PLAN
    vector<int> node_order = {2,1,0};
      auto best_plan_tree = new PlanTree(node_order);
#else
    auto best_plan_tree = (new PlanGenerator(kv_store_, bgp_query.get(), statistics, var_candidates_cache))->get_random_plan();
    // get_plan(basic_query_pointer, this->kv_store_, var_candidates_cache);
#endif
    long t4 = Util::get_cur_time();

    auto bfs_result = this->ExecutionBreathFirst(bgp_query,query_info,best_plan_tree->root_node,var_candidates_cache);


    auto bfs_table = get<1>(bfs_result);
    auto pos_var_mapping = bfs_table.pos_id_map;
    auto var_pos_mapping = bfs_table.id_pos_map;
    r = std::move(bfs_table);
  }

  decltype(query_info.limit_num_) deleted_num = 0;
  auto &result_list = r.values_;

  auto var_num = r.pos_id_map->size();
  vector<double> scores_list_;
  RankTheResult(var_coefficients,bgp_query,r, scores_list_,query_info.limit_num_);

  auto t1 = Util::get_cur_time();
  UTIL_current_query_time = t1-t0;


#ifdef TOPK_DEBUG_RESULT_INFO
  vector<bool> is_predicate_var(var_num,false);
  for(unsigned i =0;i<var_num;i++)
    is_predicate_var[i] = bgp_query->get_vardescrip_by_index(i)->var_type_ == VarDescriptor::VarType::Predicate;

  auto it = result_list->begin();
  for(decltype(result_list->size()) i =0;i<result_list->size();i++)
  {
    auto rec = *it;
    cout<<" record["<<i<<"]"<<" score:"<<scores_list_[i];
    test_scores.push_back(scores_list_[i]);
    result_ids.push_back(*rec);
    for(unsigned j =0;j<var_num;j++)
      if(is_predicate_var[(*r.pos_id_map)[j]])
        cout<<" "<<kv_store_->getPredicateByID((*rec)[j]);
      else
        cout<<" "<<kv_store_->getStringByID((*rec)[j]);
    cout<<endl;
    it++;
  }

#endif
  return std::make_tuple(true,r);
}


void Optimizer::RankTheResult(std::shared_ptr<std::unordered_map<std::string,double>> &var_coefficient, shared_ptr<BGPQuery> bgp_query,
                              IntermediateResult &result, vector<double>& score_list,int k ) {

  auto var_num = result.pos_id_map->size();
  vector<double> pos_coefficient(var_num);
  for(decltype(var_num) pos = 0; pos< var_num;pos++)
  {
    pos_coefficient[pos] = 0;
    auto var_id = (*result.pos_id_map)[pos];
    auto coefficient_it = var_coefficient->find(bgp_query->get_var_name_by_id((var_id)));
    bool has_coefficient = coefficient_it != var_coefficient->end();
    double coefficient = has_coefficient ? (*coefficient_it).second : 0.0;
    pos_coefficient[pos] = coefficient;
  }
  auto get_pos_node_score = [&pos_coefficient , this](int pos, unsigned literal_id){
    auto coefficient = pos_coefficient[pos];
    if(coefficient == 0)
      return make_pair(true,0.0);
    auto literal_string = this->kv_store_->getLiteralByID(literal_id);
    pair<bool, double> check_result =  Util::checkGetNumericLiteral(literal_string);
    bool delete_it = !get<0>(check_result);
    auto numeric_score = get<1>(check_result);
    if(!delete_it)
      return  make_pair(true,coefficient *numeric_score);
    else
      return make_pair(false,0.0);
  };

  auto old_result_list = result.values_;
  auto result_size = old_result_list->size();
  vector< pair<double,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>>> >  rank_table;
  rank_table.reserve(result_size);
  score_list.clear();
  score_list.reserve(k);
  for(auto it = old_result_list->begin(); it != old_result_list->end(); it = old_result_list->erase(it))
  {
    double tuple_score = 0;
    bool flag = true;
    for(decltype(var_num) pos = 0; pos< var_num;pos++)
    {
      auto node_score_pair = get_pos_node_score(pos,(**it)[pos]);
      if(!node_score_pair.first)
        flag = false;
      tuple_score += node_score_pair.second;
    }
    if(flag)
      rank_table.push_back(make_pair(tuple_score, std::move(*it)));
  }
  old_result_list->clear();
  sort(rank_table.begin(),rank_table.end(),[](const pair<double,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>>>& a,
  const pair<double,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>>>& b){
    return a.first < b.first;
  });

  auto &new_table = result.values_;
  for(int i = 0;i<rank_table.size() && i<k ;i++) {
    new_table->push_back(std::move(rank_table[i].second));
    score_list.push_back(rank_table[i].first);
  }
}

#endif