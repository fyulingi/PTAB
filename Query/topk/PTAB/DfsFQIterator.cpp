
#include "DfsFQIterator.h"
#include "DfsFRIterator.h"
#include "DfsOWIterator.h"

void DfsFQIterator::TryGetNext(DfsHelpInfo *dfs_help_info)
{
  // The preparation of GetFirst
  if (this->not_fully_explored_) {
    double cost = 0;
    for(auto& fr_ow_it:this->fr_ow_iterators_) {
      fr_ow_it->TryGetNext(dfs_help_info);
      if(fr_ow_it->pool_.empty())
        return;
       cost += fr_ow_it->pool_[0].cost;
    }
    // assemble the first element
    auto e_seq = sequence(fr_ow_iterators_.size(), 0);
    // this->dynamic_trie_.insert(e_seq);
    this->compressed_vector_.Insert(e_seq);
    this->queue_.push(DfsFqElement(std::move(e_seq),DfsValue::Exact,cost));
    this->not_fully_explored_ = false;
  }

  while (!this->queue_.empty()) {
    auto top_element = this->queue_.popMin();
    auto type_value = top_element.value_type;
    if (type_value == DfsValue::Range) {
      bool success_change = true;
      double cost = 0.0;
      // to make it a Exact match
      for (unsigned int i = 0;i<this->fr_ow_iterators_.size();i++) {
        auto &fr_ow_it = this->fr_ow_iterators_[i];
        auto seq_id = top_element.seq[i];
        // means the uncertainty comes from this branch
        if (fr_ow_it->pool_.size() <= seq_id)
          fr_ow_it->TryGetNext(dfs_help_info);
        // if not success
        if (fr_ow_it->pool_.size() <= seq_id) {
          success_change = false;
          break;
        }
        else
          cost += fr_ow_it->pool_[seq_id].cost;
      }
//#ifdef ITERATOR_COUNT
//      search_node_count += fr_ow_iterators_.size();
//#endif
      if(!success_change)
        continue;
      top_element.cost = cost;
      top_element.value_type = DfsValue::Exact;
      this->queue_.push(std::move(top_element));
    } else if(type_value == DfsValue::Exact) {
      auto cost = top_element.cost;
      this->e_pool_.push_back(FqElement(top_element.seq,cost));
      this->pool_.push_back(element{0,static_cast<unsigned int>(this->pool_.size()-1),cost});
      // and now, it is time to calculate it descendents
      auto &seq = top_element.seq;
      this->compressed_vector_.Insert(seq);
      for(unsigned int j=0; j<this->fr_ow_iterators_.size(); j++) {
        seq[j] += 1;
        //if(this->dynamic_trie_.detect(seq))
        if(this->compressed_vector_.AllParentsInserted(seq)) {
          auto& fr_ow_it = this->fr_ow_iterators_[j];
          if(fr_ow_it->pool_.size() == seq[j] && !fr_ow_it->Exhausted()) {
            auto estimated_cost =  top_element.cost - fr_ow_it->pool_[seq[j] - 1].cost;
            // OWs don't have estimation
            if(fr_ow_it->Type() == OrderedListType::FR)
            {
              fr_ow_it->Estimate();
              estimated_cost += fr_ow_it->estimate_[0];
              this->queue_.push(DfsFqElement(seq,DfsValue::Range,estimated_cost));
            }
          }
          else if(fr_ow_it->pool_.size() > seq[j]) {
            auto cost = top_element.cost - fr_ow_it->pool_[seq[j] - 1].cost + fr_ow_it->pool_[seq[j]].cost;
            this->queue_.push(DfsFqElement(seq,DfsValue::Exact,cost));
          }
        }
        seq[j] -= 1;
      }
      // update the estimate
      if(this->queue_.empty())
        this->SetExhausted(true);
      else
        this->estimate_[0] = this->queue_.findMin().cost;
      // Terminate Once Found
      return;
    }
    else
      throw string("wrong DfsValue in DfsFQIterator::TryGetNext");
  }
}

void DfsFQIterator::Insert(std::shared_ptr<DfsList> FR_OW_iterator) {
  this->fr_ow_iterators_.push_back(std::move(FR_OW_iterator));
}

/**
 * Insert a bulk of FR or OW iterators.
 * inserting one certain type each time each time
 * certain type [i] specified by it's the i-th child of its father
 * @param FR_OW_iterators
 */
void DfsFQIterator::Insert(std::vector<std::shared_ptr<DfsList>> FR_OW_iterators) {
  this->fr_ow_iterators_ = std::move(FR_OW_iterators);
}

void DfsFQIterator::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                              NodeOneChildVarPredicatesPtr predicate_information) {
  record->push_back(this->node_id_);
  auto &seq = this->e_pool_[i_th].seq;
  for(unsigned int i =0; i<this->fr_ow_iterators_.size(); i++) {
    if(fr_ow_iterators_[i]->Type() ==OrderedListType::OW) {
      auto ow_predicates = this->types_predicates_[i];
      fr_ow_iterators_[i]->GetResult(seq[i], record,ow_predicates);
    }
    else
      fr_ow_iterators_[i]->GetResult(seq[i], record);
  }
}

void DfsFQIterator::Estimate() {
  // only calculate the estimate_[0] with descendents' estimate_[0]
  if (!this->not_fully_explored_) {
    this->estimate_[0] = 0;
    for (const auto &descent_it : this->fr_ow_iterators_)
      this->estimate_[0] += descent_it->estimate_[0];
    return;
  }

  // initial
  for (unsigned int i = 0; i < D_HOP_LOG; i++)
    this->estimate_[i] = 0;

  for (const auto &descent_it : this->fr_ow_iterators_) {
    if (descent_it->Type() == OrderedListType::OW)
      continue;
    for (unsigned int i = 0; i < D_HOP_LOG - 1; i++) {
      if (descent_it->pool_.empty())
        this->estimate_[i] += descent_it->estimate_[i + 1];
      else
        this->estimate_[i] += std::min(descent_it->estimate_[i + 1], descent_it->pool_[0].cost);
    }
  }
  for (const auto &descent_it : this->fr_ow_iterators_) {
    if (descent_it->Type() != OrderedListType::OW)
      continue;
    for (unsigned int i = 0; i < D_HOP_LOG; i++) {
      if (descent_it->pool_.empty())
        this->estimate_[i] += descent_it->estimate_[i];
      else
        this->estimate_[i] += std::min(descent_it->estimate_[i], descent_it->pool_[0].cost);
    }
  }
}
