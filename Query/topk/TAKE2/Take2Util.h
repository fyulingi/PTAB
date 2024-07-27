
#ifndef TAKE2_UTIL_H_
#define TAKE2_UTIL_H_


#include "../TopKUtil.h"
#include "./Take2Subspace.h"

namespace Take2Util{

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<Take2FQIterator>>  AssemblingFrOw(set<TYPE_ENTITY_LITERAL_ID> &fq_ids,
                                                                                std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>>& node_scores, int k,
                                                                                vector<std::map<TYPE_ENTITY_LITERAL_ID, std::shared_ptr<Take2FRIterator>>> &descendents_FRs,
                                                                                std::vector<std::map<TYPE_ENTITY_LITERAL_ID,std::pair<std::shared_ptr<Take2OWIterator>, NodeOneChildVarPredicatesPtr>>>&descendents_OWs);


std::map<TYPE_ENTITY_LITERAL_ID, // parent id
         std::pair<std::shared_ptr<Take2OWIterator>, // its OW
                   NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
GenerateOWs(int child_var,std::shared_ptr<TopKPlanUtil::TreeEdge> &tree_edges_,
            std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
            std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
            TopKUtil::Env *env);

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<Take2FRIterator>>  GenerateFRs(int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge>& tree_edges_,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                                                                             TopKTreeNode *child_tree_node, TopKUtil::Env *env);

std::shared_ptr<Take2FRIterator> BuildIteratorTree(const shared_ptr<TopKSearchPlan> tree_search_plan, TopKUtil::Env *env);

std::unique_ptr<Take2Subspace> WrapStartSpace(std::shared_ptr<Take2FRIterator> &root_fr);
std::shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>> Divide(unsigned int k,
                                                       const std::unique_ptr<Take2Subspace>& s_ptr,
                                                       Take2SpaceHeap *collection);
}
#endif //TAKE2_UTIL_H_
