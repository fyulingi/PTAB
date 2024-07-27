
#include "DfsFQIterator.h"
#include "DfsFRIterator.h"
#include "DfsOWIterator.h"
#include "DfsUtil.h"
#include "DfsUtilDynamic.h"
/**
 * It's a complicated task
 * Using heap to prune twig
 * @param k
 */
void DfsFRIterator::TryGetNext(DfsHelpInfo* dfs_help_info)
{
  if(this->not_fully_explored_)
  {
    this->ExploreAllDescendents(dfs_help_info);
    this->not_fully_explored_ = false;
  }

  // fully explored
  while(!queue_.empty())
  {
    auto min_e = this->queue_.popMin();
    auto node_id = min_e.node;
    auto index = min_e.index;
    // auto &fq = this->fqs_map_[node_id];
    auto &fq = min_e.fq_;
    auto &predicates = min_e.predicates_;
    // If we already find one exact match, then we happily return the result
    if(min_e.value_type==DfsValue::Exact)
    {
      this->pool_.push_back(element{.node = node_id, .index = index, .cost = min_e.cost});
      this->fqs_.push_back(min_e.fq_);
      this->predicates_vec.push_back(min_e.predicates_);
      auto next_index = index + 1;
      // If the fq
      if(fq->pool_.size() > next_index)
      {
        this->queue_.push(DfsFrElement{.node=node_id,
            .index= next_index,
            .value_type=DfsValue::Exact,
            .cost = fq->pool_[next_index].cost,
            .fq_ = fq,
            .predicates_ = predicates});
      }
      else if(fq->pool_.size() == next_index &&  !fq->Exhausted())
      {
        this->queue_.push(DfsFrElement{.node=node_id,
            .index=next_index,
            .value_type=DfsValue::Range,
            .cost = fq->estimate_[0],
            .fq_ = fq,
            .predicates_ = predicates});
      }
      if(this->queue_.empty())
        this->SetExhausted(true);
      return;
    }
    // If it is a range match, try to explore it
    if (fq->pool_.size() < index + 1)
      fq->TryGetNext(dfs_help_info);
    // successfully find one, add to the queue
    if(fq->pool_.size() >= index + 1)
    {
      this->queue_.push(DfsFrElement{.node=node_id,
          .index=index,
          .value_type=DfsValue::Exact,
          .cost = fq->pool_[index].cost,
          .fq_ = fq,
          .predicates_ = predicates});
    }
  }

  if(this->queue_.empty())
    this->SetExhausted(true);
  else// because the pool size increase, the 'estimate' should be updated
    this->estimate_[0] = this->queue_.findMin().cost;
  // If we explore all element in queue, then return
  return;

}

/**
 * Insert One FQ iterator into the FR
 * @param k
 * @param fq_pointer
 */
void DfsFRIterator::Insert(TYPE_ENTITY_LITERAL_ID fq_id,
                           std::shared_ptr<DfsList> fq_pointer,
                           OnePointPredicatePtr predicates_vec) {
  // (*this->type_predicates_)[fq_id] = std::move(predicates_vec);
  DfsFrElement e{.node=fq_id, .index=0};
  if (fq_pointer->pool_.empty()) {
    e.cost = fq_pointer->estimate_[0];
    e.value_type = DfsValue::Range;
  } else {
    e.cost = fq_pointer->pool_[0].cost;
    e.value_type = DfsValue::Exact;
  }
  e.fq_ = std::move(fq_pointer);
  e.predicates_ = std::move(predicates_vec);
  // this->fqs_map_[fq_id] = std::move(fq_pointer);
  queue_.push(e);
}

/**
 *
 * @param i_th
 * @param record
 * @param predicate_information not used, because FR itself saves the information
 */
void DfsFRIterator::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                              NodeOneChildVarPredicatesPtr predicate_information) {
  auto fq = this->pool_[i_th];
  auto fq_i_th = fq.index;
  auto fq_id = fq.node;
  auto fq_pointer = this->fqs_[i_th];
  // auto fq_pointer = this->fqs_map_[fq_id];
  if(this->predicates_vec[i_th] != nullptr)
    for (auto pre_id:*this->predicates_vec[i_th])
      record->push_back(pre_id);
//  if(this->type_predicates_->find(fq_id)!=this->type_predicates_->end()) {
//    auto predicates = (*this->type_predicates_)[fq_id];
//    for (auto pre_id:*predicates)
//      record->push_back(pre_id);
//  }
  fq_pointer->GetResult(fq_i_th,record);
}


DfsFRIterator::DfsFRIterator() {
  // this->type_predicates_=std::make_shared<NodeOneChildVarPredicates>();
}

void DfsFRIterator::Estimate() {
  // the 'this->estimate_' vector has been calculated already
  if(this->has_estimated_)
    return;
  const auto & min_e = this->queue_.findMin();
  // const auto& fq = this->fqs_map_[min_e.node];
  auto &min_fq = min_e.fq_;
  min_fq->Estimate();
  std::copy_n(min_fq->estimate_, D_HOP_LOG, this->estimate_);
}


void DfsFRIterator::ExploreAllDescendents(DfsHelpInfo* dfs_help_info) {
  int fq_var_id;
  //fq_var_id = this->fqs_map_.begin()->second->GetVarId();
  if(!this->queue_.empty())
  {
    fq_var_id = this->queue_.findMin().fq_->GetVarId();
  }
  else if(!this->fqs_.empty())
    fq_var_id = this->fqs_[0]->GetVarId();
  else
    return;

  auto fq_it = this->unexplored_fqs_->Begin();
  while(fq_it !=  this->unexplored_fqs_->End()) {
    auto chosen_id_to_explore = fq_it->first;
    size_t useless_distinct_level = 0;
    size_t useless_uncertain_level = 0;
    if(dfs_help_info->DfsNew)
    {
      auto fq = DfsUtilCompressedVector::ExploreFQ(chosen_id_to_explore,
                                                   fq_var_id,
                                                   useless_distinct_level,
                                                   useless_uncertain_level,
                                                   dfs_help_info);
      if (fq != nullptr)
        this->Insert(chosen_id_to_explore, fq, fq_it->second);

    }
    else {
      auto fq = DfsUtilDynamic::ExploreFQ(chosen_id_to_explore,
                                          fq_var_id,
                                          useless_distinct_level,
                                          useless_uncertain_level,
                                          dfs_help_info);
      if (fq != nullptr)
        this->Insert(chosen_id_to_explore, fq, fq_it->second);
    }
    fq_it = this->unexplored_fqs_->Erase(fq_it);
  }
  this->Estimate();
}
