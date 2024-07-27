/*=============================================================================
# Filename: Optimizer.h
# Author: Lei Yang, Yuqi Zhou
# Mail: yangop8@pku.edu.cn. zhouyuqi@pku.edu.cn
=============================================================================*/

#ifndef _DATABASE_OPTIMIZER_H
#define _DATABASE_OPTIMIZER_H


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
#include "../Util/OrderedVector.h"
#include "./Executor.h"
#ifdef TOPK_SUPPORT
#include "../Query/topk/DPB/DPBTopKUtil.h"
#include "../Query/topk/DPP/DPPTopKUtil.h"
#include "../Query/topk/kTPM/kTPMUtil.h"
#include "../Query/topk/SAE/SAEUtil.h"
#include "../Query/topk/PTAB/DfsUtil.h"
#include "../Query/topk/PTAB/DfsUtilDynamic.h"
#include "../Query/topk/TAKEALL/TakeAllUtil.h"
#include "../Query/topk/TAKE2/Take2Util.h"
#include "../Query/topk/EAGER/EagerUtil.h"
#include "../Query/topk/PTAB/CompressedVector.h"

#endif // TOPK_SUPPORT
#include <unordered_map>
#include <map>
#include <cstring>
#include <climits>
#include <algorithm>
#include <list>
#include <cstdlib>
#include <cmath>
#include <random>

//#define FEED_PLAN
#define OPTIMIZER_DEBUG_INFO
using namespace std;

enum class TopKStrategy{
  DP_B,
  DP_P,
  K_TPM,
  SAE,
  DFS_OLD,
  PTAB,
  TAKE_ALL,
  TAKE2,
  EAGER,
  RankAfterMatching
};

using TableContent = list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>;
using TableContentShardPtr = shared_ptr<list<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>;
using PositionValue = std::map<TYPE_ENTITY_LITERAL_ID, TYPE_ENTITY_LITERAL_ID>;
using PositionValueSharedPtr = std::shared_ptr<PositionValue>;
using QueryPlanSharedPtr = shared_ptr<QueryPlan>;

enum class BasicQueryStrategy{
  Normal, // method 0
  Special,// method 1-5
  limitK,
  TopK
};

struct QueryInfo{
  bool limit_;
  int limit_num_;
  bool is_distinct_;
  shared_ptr<vector<QueryTree::Order>> ordered_by_expressions_;
};
/*
 *
 * 还差一个 copyToResult 的操作
 *
1. Optimizer类在generalevaluation里面初始化，并且在generalevaluation调用do_query()
----------
2. do_query()首先判断handler0或者handler1-5
3. 在handler0里面首先调用enum_query_plan()，按照BFS或者DFS生成所有执行计划
4. 生成完执行计划后调用choose_exec_plan()选择一个计划执行
5. 在choose_exec_plan()中每一个执行计划需要层次调用cost_model()->cardinality_estimator()->update_cardinality_cache()
 ----------
6. 按照点的顺序首先需要对起始点调用gen_filter获得该点的候选集 [done]
    而后的每一个点调用join_one_step()
7. 在每一步的join_one_step()中需要根据queryplan调用三种基本join之一，同时需要update中间结果，并且每一步的开始需要check一下中间结果能不能利用
8. 每一步都做完之后，copytoResult()
*/

class Optimizer
{
 public:

  Optimizer(KVstore* kv_store, Statistics *statistics, TYPE_TRIPLE_NUM* pre2num, TYPE_TRIPLE_NUM* pre2sub,
             TYPE_TRIPLE_NUM* pre2obj, TYPE_PREDICATE_ID limitID_predicate, TYPE_ENTITY_LITERAL_ID limitID_literal,
             TYPE_ENTITY_LITERAL_ID limitID_entity,  HopIndex* hop_index, shared_ptr<Transaction> txn);
  ~Optimizer()=default;

  tuple<bool,shared_ptr<IntermediateResult>> DoQuery(SPARQLquery&,QueryInfo); // the whole process
  tuple<bool,shared_ptr<IntermediateResult>> DoQuery(std::shared_ptr<BGPQuery> bgp_query,QueryInfo query_info,
                                                     TopKStrategy top_k_method=TopKStrategy::DP_B); // the whole process

 private:
  tuple<bool,shared_ptr<IntermediateResult>> MergeBasicQuery(SPARQLquery &sparql_query);

    /* the function is not implemented yet
   TODO: re-choose plan in every iteration
    implement a method to feed plan manually first
    tuple<bool,shared_ptr<IntermediateResult>> execution(BasicQuery*, shared_ptr<QueryPlan>);*/


  BasicQueryStrategy ChooseStrategy(BasicQuery *basic_query,QueryInfo *query_info);
  BasicQueryStrategy ChooseStrategy(std::shared_ptr<BGPQuery> bgp_query,QueryInfo *query_info);

  std::shared_ptr<std::vector<IntermediateResult>> GenerateResultTemplate(shared_ptr<QueryPlan> query_plan);

  tuple<bool,IntermediateResult> ExecutionDepthFirst(shared_ptr<BGPQuery>  bgp_query,
                                                       shared_ptr<QueryPlan> query_plan,
                                                       QueryInfo query_info,
                                                       PositionValueSharedPtr id_pos_mapping);


  tuple<bool,IntermediateResult> DepthSearchOneLayer(shared_ptr<QueryPlan> query_plan,
                                                     int layer_count,
                                                     int &result_number_till_now,
                                                     int limit_number,
                                                     TableContentShardPtr table_content_ptr,
                                                     IDCachesSharePtr id_caches,
                                                     shared_ptr<vector<IntermediateResult>> table_template);



#ifdef TOPK_SUPPORT

  tuple<bool,IntermediateResult> ExecutionTopK_B(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                 const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                 std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_P(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                 const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                 std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_kTPM(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                    const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                    std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_SAE(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                    const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                    std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_DfsOld(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                      const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                      std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_Dfs(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                   const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                   std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_TakeAll(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                       const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                       std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_Take2(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                     const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                     std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_Eager(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                     const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                     std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

  tuple<bool,IntermediateResult> ExecutionTopK_RAM(shared_ptr<BGPQuery> bgp_query, shared_ptr<TopKSearchPlan> &tree_search_plan,
                                                     const QueryInfo& query_info, IDCachesSharePtr id_caches,
                                                     std::shared_ptr<std::unordered_map<std::string,double>> var_coefficients);

#endif
  /*copy the result to vector<unsigned*> & */
  bool CopyToResult(vector<unsigned*> *target, BasicQuery *basic_query, const shared_ptr<IntermediateResult>& result);

  bool CopyToResult(vector<unsigned*> *target, shared_ptr<BGPQuery> bgp_query, IntermediateResult result);

  void RankTheResult(std::shared_ptr<std::unordered_map<std::string,double>> &var_coefficient, shared_ptr<BGPQuery> bgp_query,
                     IntermediateResult &result, vector<double>& score_list,int k);

  void Cartesian(int, int,int,unsigned*,const shared_ptr<vector<Satellite>>&,vector<unsigned*>*,BasicQuery *);

  tuple<bool,PositionValueSharedPtr, TableContentShardPtr> ExecutionBreathFirst(BasicQuery* basic_query,const QueryInfo& query_info,Tree_node* plan_tree,IDCachesSharePtr id_caches);

  tuple<bool,IntermediateResult> ExecutionBreathFirst(shared_ptr<BGPQuery> bgp_query,
                                                      QueryInfo query_info,
                                                      Tree_node* plan_tree_node,
                                                      IDCachesSharePtr id_caches);
 private:
  KVstore* kv_store_;
  Statistics* statistics;
  HopIndex* hop_index_;
  shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>> order_by_list_; // empty if not using 'orderby'
  TYPE_ENTITY_LITERAL_ID limit_num_; // -1 if not limit result size

  bool is_edge_case; //Strategy 1-6 or Strategy 0
  int current_basic_query_; // updated by result_list.size()
  shared_ptr<vector<BasicQuery*>> basic_query_list_; //fork from SPARQLQuery, I dont know why
  shared_ptr<vector<tuple<BasicQuery*, shared_ptr<vector<QueryPlan>>>>> candidate_plans_;
  shared_ptr<vector<QueryPlan>> execution_plan_;
  shared_ptr<vector<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>> result_list_; // vector<unsigned*>* result_list;

  //TODO: shared_ptr<BasicQuery> may brings wrong matching
  shared_ptr<vector<map<BasicQuery*,vector<shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>>> join_cache_; // map(sub-structure, result_list)
  shared_ptr<vector<map<BasicQuery*,shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>>>>> cardinality_cache_; // map(sub-structure, cardinality), not in statistics

  FILE* fp_;
  TYPE_TRIPLE_NUM* pre2num_;
  TYPE_TRIPLE_NUM* pre2sub_;
  TYPE_TRIPLE_NUM* pre2obj_;
  TYPE_PREDICATE_ID limitID_predicate_;
  TYPE_ENTITY_LITERAL_ID limitID_literal_;
  TYPE_ENTITY_LITERAL_ID limitID_entity_;
  shared_ptr<Transaction> txn_;

  //静态情况下可以存储weight
  shared_ptr<vector<double>> weight_list_;
  Executor executor_;

  tuple<bool,bool,bool> PrepareInitial(shared_ptr<BGPQuery> bgp_query,
                            shared_ptr<FeedOneNode> join_a_node_plan) const;
};


#endif