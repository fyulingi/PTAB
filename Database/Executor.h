/*=============================================================================
# Filename: Executor.h
# Author: Yuqi Zhou
# Mail: zhouyuqi@pku.edu.cn
# Last Modified:  2021/8/10.
=============================================================================*/
#ifndef GSTOREGDB_DATABASE_EXECUTOR_H_
#define GSTOREGDB_DATABASE_EXECUTOR_H_
#include "../Util/Util.h"
#include "../Query/SPARQLquery.h"
#include "../Query/BasicQuery.h"
#include "../Query/IDList.h"
#include "../KVstore/KVstore.h"
#include "../VSTree/VSTree.h"
#include "TableOperator.h"
#include "ResultTrigger.h"
#include "../Query/QueryPlan.h"
#include "../Query/QueryTree.h"
#include "Join.h"
#include "./Statistics.h"
#include "./PlanTree.h"
#include "./PlanGenerator.h"

using namespace std;

template <typename Container> // we can make this generic for any container [1]
struct container_hash{
  std::size_t operator()(Container const& c) const {
    return boost::hash_range(c.begin(), c.end());
  }
};

/**
 * Operation for Table is doing here.
 * any operation to write id_pos/pos_id map should be in here.
 * IntermediateResult is the standard format Executor communicates
 * with others
 */
class Executor {
 private:
  KVstore *kv_store_;
  std::shared_ptr<Transaction> txn_;
  TYPE_PREDICATE_ID limitID_predicate_;
  TYPE_ENTITY_LITERAL_ID limitID_literal_;
  TYPE_ENTITY_LITERAL_ID limitID_entity_;


  tuple<bool,TableContentShardPtr> OneEdgeConstraintFilter(EdgeInfo edge_info,
                                                           EdgeConstantInfo edge_table_info,
                                                           TableContentShardPtr table_content_ptr,
                                                           PositionValueSharedPtr id_pos_mapping,
                                                           IDCachesSharePtr id_caches);

  tuple<bool,TableContentShardPtr> FilterAVariableOnIDList(shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>> candidate_list,
                                                           TYPE_ENTITY_LITERAL_ID var_id,
                                                           TableContentShardPtr table_content_ptr,
                                                           PositionValueSharedPtr id_posing_mapping);

  TableContentShardPtr ConvertToTable(std::shared_ptr<IDList> id_list);

  bool AddConstantCandidates(EdgeInfo edge_info,TYPE_ENTITY_LITERAL_ID targetID,IDCachesSharePtr id_caches);

  std::tuple<bool, TableContentShardPtr> GetAllSubObjId(bool need_literal);
  std::tuple<bool, TableContentShardPtr> GetAllPreId();
  std::shared_ptr<IDList> ExtendRecordOneNode(shared_ptr<FeedOneNode> one_step_join_node_,
                                              PositionValueSharedPtr id_pos_mapping,
                                              IDCachesSharePtr id_caches,
                                              TYPE_ENTITY_LITERAL_ID new_id,
                                              list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>::iterator record_iterator) const;

  std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> ExtendRecordTwoNode(shared_ptr<FeedTwoNode> one_step_join_node_,
                                                                           PositionValueSharedPtr id_pos_mapping,
                                                                           IDCachesSharePtr id_caches,
                                                                           TYPE_ENTITY_LITERAL_ID new_id1,
                                                                           TYPE_ENTITY_LITERAL_ID new_id2,
                                                                           list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>::iterator record_iterator) const;

  shared_ptr<IDList> CandidatesWithConstantEdge(shared_ptr<vector<EdgeInfo>> edge_info_vector) const;

 public:
  Executor(KVstore *kv_store,std::shared_ptr<Transaction> txn,TYPE_PREDICATE_ID limitID_predicate,
           TYPE_ENTITY_LITERAL_ID limitID_literal,TYPE_ENTITY_LITERAL_ID limitID_entity):
           kv_store_(kv_store),txn_(txn),limitID_predicate_(limitID_predicate),limitID_literal_(limitID_literal),
           limitID_entity_(limitID_entity){};

  tuple<bool, IntermediateResult> JoinANode(IntermediateResult old_table,
                                            IDCachesSharePtr id_caches,
                                            shared_ptr<FeedOneNode> feed_plan);

  std::tuple<bool,IntermediateResult> InitialTableOneNode(std::shared_ptr<FeedOneNode> feed_plan,
                                                          bool is_entity,
                                                          bool is_literal,
                                                          bool is_predicate,
                                                          IDCachesSharePtr id_caches);

  std::tuple<bool,IntermediateResult> InitialTableTwoNode(std::shared_ptr<FeedTwoNode> join_plan,
                                                          IDCachesSharePtr id_caches);

  std::tuple<bool, IntermediateResult> JoinTwoNode(shared_ptr<FeedTwoNode> join_two_node_,
                                                   IntermediateResult old_table,
                                                   IDCachesSharePtr id_caches);

  tuple<bool,IntermediateResult> JoinTable(shared_ptr<JoinTwoTable> join_plan,
                                           IntermediateResult table_a,
                                           IntermediateResult table_b);

  tuple<bool,IntermediateResult> ANodeEdgesConstraintFilter(shared_ptr<FeedOneNode> check_plan,
                                                            IntermediateResult old_table,
                                                            IDCachesSharePtr id_caches);

  bool CacheConstantCandidates(shared_ptr<FeedOneNode> one_step, IDCachesSharePtr id_caches);

  static void UpdateIDList(const std::shared_ptr<IDList>& valid_id_list, unsigned* id_list, unsigned id_list_len,bool id_list_prepared);

  bool UpdateIDCache(std::shared_ptr<StepOperation> step_operation, IDCachesSharePtr id_caches);


};

#endif //GSTOREGDB_DATABASE_EXECUTOR_H_
