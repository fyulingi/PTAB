
#ifndef KTPM_UTIL_H_
#define KTPM_UTIL_H_


#include "../TopKUtil.h"
#include "./Subspace.h"

namespace kTPMUtil{

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<kTPMFQIterator>>  AssemblingFrOw(set<TYPE_ENTITY_LITERAL_ID> &fq_ids,
                                                                                std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>>& node_scores, int k,
                                                                                vector<std::map<TYPE_ENTITY_LITERAL_ID, std::shared_ptr<kTPMFRIterator>>> &descendents_FRs,
                                                                                std::vector<std::map<TYPE_ENTITY_LITERAL_ID,std::pair<std::shared_ptr<kTPMOWIterator>, NodeOneChildVarPredicatesPtr>>>&descendents_OWs);


std::map<TYPE_ENTITY_LITERAL_ID, // parent id
         std::pair<std::shared_ptr<kTPMOWIterator>, // its OW
                   NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
GenerateOWs(int child_var,std::shared_ptr<TopKPlanUtil::TreeEdge> &tree_edges_,
            std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
            std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
            TopKUtil::Env *env);

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<kTPMFRIterator>>  GenerateFRs(int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge>& tree_edges_,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                                                                             TopKTreeNode *child_tree_node, TopKUtil::Env *env);

std::shared_ptr<kTPMFRIterator> BuildIteratorTree(const shared_ptr<TopKSearchPlan> tree_search_plan, TopKUtil::Env *env);

std::unique_ptr<Subspace> WrapStartSpace(std::shared_ptr<kTPMFRIterator> &root_fr);
std::shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>> Divide(unsigned int k,
                                                       const std::unique_ptr<Subspace>& s_ptr,
                                                       SpaceHeap *collection);
}
#endif //KTPM_UTIL_H_
