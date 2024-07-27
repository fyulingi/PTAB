
#include "DPBanyList.h"

// getFirst() also implemented here
void DPBanyFR::TryGetNext(unsigned int k) {
  // get first
  if(this->pool_.size()==0)
  {
    auto queue_min = this->queue_.findMin();
    element e;
    e.cost = queue_min.cost;
    e.index = queue_min.index;
    e.node = queue_min.node;
    this->pool_.push_back(e);
    this->pool_fqs_.push_back(queue_min.fq_);
    this->pool_predicates_.push_back(queue_min.predicates_);
    return;
  }
  if(this->queue_.empty())
    return;
  auto m = this->pool_.size();
  auto em = this->pool_.back();
  // auto em_pointer = this->fqs_map_[em.node];
  auto em_pointer = this->pool_fqs_.back();
  auto em_predicates = this->pool_predicates_.back();
  this->queue_.popMin();
  if(NextEPoolElement(k,em_pointer,em.index+1))
  {
    DpbAnyElement e;
    e.node = em.node;
    e.index = em.index + 1;
    e.cost = em.cost + this->DeltaCost(em_pointer,e.index);
    e.fq_ = em_pointer;
    e.predicates_ = em_predicates;
    this->queue_.push(e);
  }
  if(!queue_.empty()) {
    if (this->queue_.size() > k - m)
      this->queue_.popMax();
    auto queue_min = this->queue_.findMin();
    element e;
    e.cost = queue_min.cost;
    e.index = queue_min.index;
    e.node = queue_min.node;
    this->pool_.push_back(e);
    this->pool_fqs_.push_back(queue_min.fq_);
    this->pool_predicates_.push_back(queue_min.predicates_);
  }
}


/**
 * Insert One FQ iterator into the FR
 * @param k
 * @param fq_pointer
 */
void DPBanyFR::Insert(unsigned int k,
                           TYPE_ENTITY_LITERAL_ID fq_id,
                           std::shared_ptr<DPBanyList> fq_pointer,
                           OnePointPredicatePtr predicates_vec) {
  // (*this->type_predicates_)[fq_id] = predicates_vec;
  auto cost = fq_pointer->pool_[0].cost;
  DpbAnyElement e{};
  e.cost = cost;
  e.index = 0;
  e.node = fq_id;
  e.fq_ = std::move(fq_pointer);
  e.predicates_ = std::move(predicates_vec);
  // this->fqs_map_[fq_id] = std::move(fq_pointer);
  queue_.push(e);
  if(queue_.size()>k)
    queue_.popMax();
}

double DPBanyFR::DeltaCost(std::shared_ptr<DPBanyList>& node_pointer, int index) {
  auto delta =  node_pointer->pool_[index].cost - node_pointer->pool_[index-1].cost;
  return delta;
}

bool DPBanyFR::NextEPoolElement(unsigned int k, std::shared_ptr<DPBanyList>& node_pointer, unsigned int index) {
  if(index == node_pointer->pool_.size())
    node_pointer->TryGetNext(k);
  if(index < node_pointer->pool_.size())
    return true;
  else
    return false;
}

/**
 *
 * @param i_th
 * @param record
 * @param predicate_information not used, because FR itself saves the information
 */
void DPBanyFR::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                           NodeOneChildVarPredicatesPtr predicate_information) {
  auto fq = this->pool_[i_th];
  auto fq_i_th = fq.index;
  auto fq_id = fq.node;
  auto fq_pointer = this->pool_fqs_[i_th];
  // this->fqs_map_[fq_id];
  // if(this->type_predicates_->find(fq_id)!=this->type_predicates_->end()) {
  if(!this->pool_predicates_.empty()) {
    auto predicates = this->pool_predicates_[i_th];
    for (auto pre_id:*predicates)
      record->push_back(pre_id);
  }
  fq_pointer->GetResult(fq_i_th,record);
}
DPBanyFR::DPBanyFR() {
  //this->type_predicates_=std::make_shared<NodeOneChildVarPredicates>();
}

void DPBanyOW::TryGetNext(unsigned int k) {
}

void DPBanyOW::Insert(unsigned int k,
    const std::vector<TYPE_ENTITY_LITERAL_ID>& ids,
                           const std::vector<double>& scores)
{
  struct ScorePair{
    TYPE_ENTITY_LITERAL_ID id;
    double cost;
    bool operator<(const ScorePair& other) const{return this->cost<other.cost;};
  };
  std::vector<ScorePair> ranks;
  this->pool_.reserve(k);
  ranks.reserve(ids.size());
  for(unsigned int  i=0;i<ids.size();i++)
    ranks.push_back(ScorePair{ids[i],scores[i]});
  std::sort(ranks.begin(),ranks.end());
  for(unsigned int i=0;i<k && i < ranks.size();i++)
  {
    element e{};
    e.index = i;
    e.cost = ranks[i].cost;
    e.node = ranks[i].id;
    this->pool_.push_back(e);
  }
}

/**
 * Default score : 0.0
 * @param ids
 */
void DPBanyOW::Insert(unsigned int k,
                           const std::vector<TYPE_ENTITY_LITERAL_ID> &ids) {
  pool_.reserve(k);
  for(unsigned int i=0;i<k && i < ids.size();i++)
  {
    element e{};
    e.index = i;
    e.cost = 0.0;
    e.node = ids[i];
    this->pool_.push_back(e);
  }
}

void DPBanyOW::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                              NodeOneChildVarPredicatesPtr predicate_information) {
  auto &i_th_element = this->pool_[i_th];
  auto node_id = i_th_element.node;
  auto predicates = (*predicate_information)[node_id];
  for(auto predicate_id:*predicates)
    record->push_back(predicate_id);
  record->push_back(node_id);
}

void DPBanyFQ::TryGetNext(unsigned int k) {
  // get first
  if(this->pool_.size()==0)
  {
    double cost = 0;
    for(unsigned int j =0; j<this->fr_ow_iterators_.size(); j++)
    {
      if(!DPBanyFQ::NextEPoolElement(k, this->fr_ow_iterators_[j], 0))
        return;
      cost += this->fr_ow_iterators_[j]->pool_[0].cost;
    }
    FqElement e;
    e.seq = sequence(fr_ow_iterators_.size(), 0);
    e.cost = cost;
    this->compressed_vector_.Insert(e.seq);
    this->queue_.push(e);

    this->e_pool_.push_back(e);
    // transfer e pool to i pool element
    element ipool_element{};
    ipool_element.cost = e.cost;
    ipool_element.index = this->pool_.size();
    //ipool_e.pointer = this;
    //this->seq_list_.push_back(e.seq);
    this->pool_.push_back(ipool_element);
    return;
  }
  if(this->queue_.empty())
    return;
  auto m = this->e_pool_.size();
  // Get Next-- insert the children of the top
  //auto m = this->e_pool_.size();
  auto em = this->e_pool_.back();
  queue_.popMin();
  auto seq = em.seq;
  for(unsigned int j=0; j<this->fr_ow_iterators_.size(); j++)
  {
    seq[j] += 1;
    if(this->compressed_vector_.AllParentsInserted(seq))
    {
      if(this->NextEPoolElement(k, this->fr_ow_iterators_[j], seq[j]))
      {
        decltype(em) ec;
        ec.cost = em.cost + this->DeltaCost(this->fr_ow_iterators_[j], seq[j]);
        ec.seq = seq;
        this->queue_.push(ec);
      }
    }
    seq[j] -= 1;
  }
  if(!queue_.empty()) {
    if (this->queue_.size() > k - m)
      this->queue_.popMax();
    auto inserted = this->queue_.findMin();
    inserted.cost += this->node_score_;
    this->e_pool_.push_back(inserted);
    // transfer e pool to i pool element
    element e{};
    e.cost = inserted.cost;
    e.index = this->pool_.size()-1;
    //this->seq_list_.push_back(std::move(inserted.seq));
    this->pool_.push_back(e);
  }
}

bool DPBanyFQ::NextEPoolElement(unsigned int k, std::shared_ptr<DPBanyList> &node_pointer, unsigned int index) {
  if(index == node_pointer->pool_.size())
    node_pointer->TryGetNext(k);
  if(index < node_pointer->pool_.size())
    return true;
  else
    return false;
}

void DPBanyFQ::Insert(std::shared_ptr<DPBanyList> FR_OW_iterator) {
  this->fr_ow_iterators_.push_back(std::move(FR_OW_iterator));
}

/**
 * Insert a bulk of FR or OW iterators.
 * inserting one certain type each time each time
 * certain type [i] specified by it's the i-th child of its father
 * @param FR_OW_iterators
 */
void DPBanyFQ::Insert(std::vector<std::shared_ptr<DPBanyList>> FR_OW_iterators) {
  this->fr_ow_iterators_ = std::move(FR_OW_iterators);
}


double DPBanyFQ::DeltaCost(std::shared_ptr<DPBanyList> &FR_OW_iterator, int index) {
  auto delta =  FR_OW_iterator->pool_[index].cost - FR_OW_iterator->pool_[index-1].cost;
  return delta;
}

void DPBanyFQ::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
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
