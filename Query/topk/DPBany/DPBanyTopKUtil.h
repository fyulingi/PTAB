
#ifndef GSTORELIMITK_QUERY_ANY_TOPK_TOPKUTIL_H_
#define GSTORELIMITK_QUERY_ANY_TOPK_TOPKUTIL_H_


#include "../TopKUtil.h"
#include "./DPBanyList.h"

namespace DPBanyUtil{

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPBanyFQ>>  AssemblingFrOw(set<TYPE_ENTITY_LITERAL_ID> &fq_ids,
                                                                                std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>> &node_scores, int k,
                                                                                vector<std::map<TYPE_ENTITY_LITERAL_ID, std::shared_ptr<DPBanyFR>>> &descendents_FRs,
                                                                                std::vector<std::map<TYPE_ENTITY_LITERAL_ID,std::pair<std::shared_ptr<DPBanyOW>, NodeOneChildVarPredicatesPtr>>>&descendents_OWs);


std::map<TYPE_ENTITY_LITERAL_ID, // parent id
         std::pair<std::shared_ptr<DPBanyOW>, // its OW
                   NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
GenerateOWs(int child_var,std::shared_ptr<TopKPlanUtil::TreeEdge>& tree_edges_,
            std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
            std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
            TopKUtil::Env *env);

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPBanyFR>>  GenerateFRs(int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge> &tree_edges_,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                                                                             TopKTreeNode *child_tree_node, TopKUtil::Env *env);

DPBanyFR* BuildIteratorTree(const shared_ptr<TopKSearchPlan> tree_search_plan, TopKUtil::Env *env);


}
#endif // GSTORELIMITK_QUERY_ANY_TOPK_TOPKUTIL_H_
