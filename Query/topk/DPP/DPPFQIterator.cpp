
#include "DPPFQIterator.h"
#include "DPPFRIterator.h"
#include "DPPOWIterator.h"

void DPPFQIterator::TryGetNext(unsigned int k,GlobalQueue *global_queue) {

  // Get Next-- insert the children of the top
  // unlike DP-B , get_first() is not necessary here
  auto em = this->e_pool_.back();
  queue_.popMin();
  auto seq = em.seq;
  for(unsigned int j=0; j<this->fr_ow_iterators_.size(); j++)
  {
    seq[j] += 1;
    if(this->dynamic_trie_.detect(seq))
    {
      if(this->NextIPoolElement(k, j,seq[j],global_queue))
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
    auto m = this->e_pool_.size();
    if (this->queue_.size() > k - m)
      this->queue_.popMax();
  }
  // corresponding to line 11 in algo 7.
  if(queue_.empty()) {
    if (this->GQ_pos_ != QueueElement::InValidPos)
      global_queue->Pop(this->GQ_pos_);
  }
  else
  {
    auto &e = this->queue_.findMin();
    auto queue_element = std::unique_ptr<QueueElement>(new QueueElement);
    queue_element->fq_iterator_ = this;
    queue_element->self_dpp_list_ = this;
    queue_element->fq_element_ =  std::unique_ptr<FqElement>(new FqElement(e));
    queue_element->type_ = QueueElement::Type::FQ;
    global_queue->Insert(std::move(queue_element));
  }
}

bool DPPFQIterator::NextEPoolElement(unsigned int k, std::shared_ptr<DPPList> node_pointer,
                                     unsigned int index,unsigned int j_th_child,GlobalQueue *global_queue) {
  if(index == node_pointer->pool_.size()) {
    if(node_pointer->Type()==OrderedListType::FR){
      auto FR_it = std::dynamic_pointer_cast<DPPFRIterator>(node_pointer);
      FR_it->TryGetNext(k,global_queue,this,j_th_child);
    }
  }
  if(index < node_pointer->pool_.size())
    return true;
  else
    return false;
}

void DPPFQIterator::Insert(std::shared_ptr<DPPList> FR_OW_iterator) {
  this->fr_ow_iterators_.push_back(FR_OW_iterator);
}

/**
 * Insert a bulk of FR or OW iterators.
 * inserting one certain type each time each time
 * certain type [i] specified by it's the i-th child of its father
 * @param FR_OW_iterators
 */
void DPPFQIterator::Insert(std::vector<std::shared_ptr<DPPList>> FR_OW_iterators) {
  this->fr_ow_iterators_ = std::move(FR_OW_iterators);
}


double DPPFQIterator::DeltaCost(std::shared_ptr<DPPList> FR_OW_iterator, int index) {
  auto delta =  FR_OW_iterator->pool_[index].cost - FR_OW_iterator->pool_[index-1].cost;
  return delta;
}

void DPPFQIterator::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
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

void DPPFQIterator::TryTriggeringRootSeq(size_t j,GlobalQueue *global_queue) {

  if(this->fr_ow_iterators_[j]->Type()==OrderedListType::OW)
  {
    auto ow_it = std::dynamic_pointer_cast<DPPOWIterator>(fr_ow_iterators_[j]);
    if(ow_it->Size()>1)
      return;// we have tried before
  }
  else if(this->fr_ow_iterators_[j]->pool_.size()>1)
    return; // we have tried before

  this->cost1_ += this->fr_ow_iterators_[j]->pool_[0].cost;
  this->branch_counter_ ++;
  if(this->branch_counter_ < this->child_type_num_) return;

  auto e = std::unique_ptr<FqElement>(new FqElement);
  e->seq = sequence(fr_ow_iterators_.size(), 0);
  e->cost = this->cost1_;
  this->dynamic_trie_.insert(e->seq);
  this->queue_.push(*e);
  GQ_Update(global_queue,move(e));
}

void DPPFQIterator::TriggerASingletonSeq(unsigned int k,size_t j_type,GlobalQueue *global_queue) {
  FqElement e;
  e.seq = sequence(fr_ow_iterators_.size(), 0);
  auto iPool_J_size = fr_ow_iterators_[j_type]->pool_.size();
  e.seq[j_type] =iPool_J_size-1;
  if(iPool_J_size == 2)
    this->sin_seq_cost_[j_type] = this->e_pool_[0].cost;
  e.cost = this->sin_seq_cost_[j_type] + DeltaCost(fr_ow_iterators_[j_type],iPool_J_size - 1);
  this->sin_seq_cost_[j_type] = e.cost;
  this->queue_.push(e);
  if(this->queue_.size() > k - this->e_pool_.size())
    this->queue_.popMax();
  auto e_unique = std::unique_ptr<FqElement>(new FqElement(e));
  GQ_Update(global_queue,std::move(e_unique));
}

bool DPPFQIterator::NextIPoolElement(unsigned int k , size_t j_th_type, unsigned int index, GlobalQueue *global_queue) {
  // line 19-25 in algo.6
  auto iPool_J = this->fr_ow_iterators_[j_th_type];
  if(iPool_J->Type() == OrderedListType::OW)
  {
    auto ow = std::dynamic_pointer_cast<DPPOWIterator>(iPool_J);
    if(index<ow->Size())  // line 19 in algo.6
      return true;
    if(index==ow->Size())
    {
      ow->TryGetNext(k);
      if(index == ow->Size())
        return false;
      auto queue_element = std::unique_ptr<QueueElement>(new QueueElement);
      queue_element->fq_iterator_ = this;
      queue_element->self_dpp_list_ = this;
      queue_element->i_element_ = std::unique_ptr<element>(new element(iPool_J->pool_[index]));
      queue_element->i_th_type = j_th_type;
      queue_element->type_ = QueueElement::Type::OW;
      global_queue->Insert(std::move(queue_element));
      return false;
    }
  }
  if(iPool_J->Type() == OrderedListType::FR)
  {
    if (index < iPool_J->pool_.size()) {
      // line 19 in algo.6
      return true;
    }
    else if(index == iPool_J->pool_.size()) {
      // line 21 in algo.6
      std::dynamic_pointer_cast<DPPFRIterator>(iPool_J)->TryGetNext(k, global_queue, this, j_th_type);
      return false;
    }
  }
  return false;
}

void DPPFQIterator::AddParent(DPPFRIterator* fq_iterator,unsigned int index)
{
  this->temp_parents.emplace_back(fq_iterator,index);
}

void DPPFQIterator::GQ_Update(GlobalQueue *global_queue, std::unique_ptr<FqElement> e) {
  // line 23-26 in algo 5.
  if( e->operator!=(this->queue_.findMin()))
    return;

  auto queue_element = std::unique_ptr<QueueElement>(new QueueElement);
  queue_element->fq_iterator_ = this;
  queue_element->self_dpp_list_ = this;
  queue_element->fq_element_= std::move(e);
  queue_element->type_ = QueueElement::Type::FQ;
  // line 25 in algo.5
  if(this->GQ_pos_ == QueueElement::InValidPos)
  {
    global_queue->Insert(std::move(queue_element));
  }
  else
  {
    global_queue->UpdateNode(std::move(queue_element),this->GQ_pos_);
  }

}
void DPPFQIterator::GetFirst(unsigned int k,GlobalQueue *global_queue) {
  double cost = 0;
  for(unsigned int j =0; j<this->fr_ow_iterators_.size(); j++)
  {
    if(!DPPFQIterator::NextEPoolElement(k, this->fr_ow_iterators_[j], 0,j,global_queue))
      return;
    cost += this->fr_ow_iterators_[j]->pool_[0].cost;
  }
  FqElement e;
  e.seq = sequence(fr_ow_iterators_.size(), 0);
  e.cost = cost;
  this->dynamic_trie_.insert(e.seq);
  this->queue_.push(e);

  this->e_pool_.push_back(e);
  // transfer e pool to i pool element
  element ipool_element{};
  ipool_element.cost = e.cost;
  ipool_element.index = this->pool_.size();
  this->pool_.push_back(ipool_element);
  return;
}
