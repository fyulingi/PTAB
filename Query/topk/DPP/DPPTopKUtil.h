
#ifndef PRIVATETOPK_QUERY_TOPK_DPP_DPPTOPKUTIL_H_
#define PRIVATETOPK_QUERY_TOPK_DPP_DPPTOPKUTIL_H_

#include "../TopKUtil.h"
#include "../CommonDP/MinMaxHeap.hpp"
#include "./DPPFQIterator.h"
#include "./DPPFRIterator.h"
#include "./DPPOWIterator.h"

namespace DPPUtil {
std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFQIterator>>  AssemblingFrOw(int fq_var,
                                                                                unsigned int no_coefficient_fr_subtree_num,
                                                                                set<TYPE_ENTITY_LITERAL_ID> &fq_ids,
                                                                                std::shared_ptr<std::map<TYPE_ENTITY_LITERAL_ID,double>> node_scores, int k,
                                                                                vector<std::map<TYPE_ENTITY_LITERAL_ID, std::shared_ptr<DPPFRIterator>>> &descendents_FRs,
                                                                                std::vector<std::map<TYPE_ENTITY_LITERAL_ID,std::pair<std::shared_ptr<DPPOWIterator>, NodeOneChildVarPredicatesPtr>>>&descendents_OWs,
                                                                                GlobalQueue *global_queue);


std::map<TYPE_ENTITY_LITERAL_ID, // parent id
         std::pair<std::shared_ptr<DPPOWIterator>, // its OW
                   NodeOneChildVarPredicatesPtr>> // predicate correspond to the OW
GenerateOWs(int parent_var,int child_var,std::shared_ptr<TopKPlanUtil::TreeEdge> tree_edges_,
            std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
            std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
            TopKUtil::Env *env);

std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFRIterator>>  GenerateFRs(int parent_var, int child_var, std::shared_ptr<TopKPlanUtil::TreeEdge> tree_edges_,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& parent_var_candidates,
                                                                             std::set<TYPE_ENTITY_LITERAL_ID>& deleted_parents,
                                                                             TopKTreeNode *child_tree_node, TopKUtil::Env *env,
                                                                             GlobalQueue *global_queue);

DPPFRIterator* BuildIteratorTree(const shared_ptr<TopKSearchPlan> tree_search_plan, TopKUtil::Env *env,GlobalQueue *global_queue);

/**
 * this function will do 2 thing:
 * 1. fill the 'parents_' field in DPPFQIterator
 * 2. insert coefficient OW into global queue
 * @param var_id if this is FQ/OW then it is the var iterator prents, else if FR, it is
 * the FQs' var id it has
 * @param dpp_list the pointer to FQ/FR/OW iterator
 * @param type FQ/FR/OW
 * @param var_tree_node the plan tree node
 * @param env
 * @param global_queue Global Queue
 * @return
 */
void RefineTree(int var_id,
                const std::shared_ptr<DPPList>& dpp_list,
                OrderedListType type,
                TopKTreeNode *var_tree_node,
                TopKUtil::Env *env,
                GlobalQueue *global_queue);

void TriggerSeq(unsigned int k,DPPFQIterator* fq_iterator,
                DPPList *fr_ow_iterator,
                unsigned int j_th_type,
                GlobalQueue *global_queue);

void RorDExpansion(unsigned int k,DPPFQIterator* fq_iterator,GlobalQueue *global_queue,TopKSearchPlan* search_plan);
};

#endif //PRIVATETOPK_QUERY_TOPK_DPP_DPPTOPKUTIL_H_
