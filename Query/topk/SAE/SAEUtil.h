
#ifndef SAE_UTIL_H_
#define SAE_UTIL_H_


#include "../TopKUtil.h"
#include "./SpaceSAE.h"

namespace SAEUtil{

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<SAEFQIterator>>  AssemblingFrOw(set<TYPE_ENTITY_LITERAL_ID> &fq_ids,
                                                                                std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>>& node_scores, int k,
                                                                                vector<std::map<TYPE_ENTITY_LITERAL_ID, std::shared_ptr<SAEFRIterator>>> &descendents_FRs,
                                                                                std::vector<std::map<TYPE_ENTITY_LITERAL_ID,std::pair<std::shared_ptr<SAEOWIterator>, NodeOneChildVarPredicatesPtr>>>&descendents_OWs);


std::map<TYPE_ENTITY_LITERAL_ID, // parent id
         std::pair<std::shared_ptr<SAEOWIterator>, // its OW
                   NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
GenerateOWs(int child_var,std::shared_ptr<TopKPlanUtil::TreeEdge> &tree_edges_,
            std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
            std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
            TopKUtil::Env *env);

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<SAEFRIterator>>  GenerateFRs(int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge>& tree_edges_,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                                                                             TopKTreeNode *child_tree_node, TopKUtil::Env *env);

std::shared_ptr<SAEFRIterator> BuildIteratorTree(const shared_ptr<TopKSearchPlan> tree_search_plan, TopKUtil::Env *env);

std::unique_ptr<SpaceSAE> WrapStartSpace(std::shared_ptr<SAEFRIterator> &root_fr);
std::shared_ptr<vector<TYPE_ENTITY_LITERAL_ID>> Divide(unsigned int k,
                                                       std::unique_ptr<SpaceSAE>& s_ptr,
                                                       SpaceHeapSAE *collection);
}
#endif //SAE_UTIL_H_
