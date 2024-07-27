
#include "GlobalQueue.h"

void GlobalQueue::HeapSwap(int i, int j)
{
  std::swap(contents_[i]->self_dpp_list_->GQ_pos_, contents_[j]->self_dpp_list_->GQ_pos_);
  std::swap(contents_[i],contents_[j]);
}

void GlobalQueue::Down(size_t pos)
{
  auto t = pos;
  auto t_value = this->contents_[t]->Score();
  auto heap_size = this->contents_.size();
  auto left_child = 2 * t;

  if(left_child>=heap_size)
    return;
  auto left_value = this->contents_[left_child]->Score();
  if(left_value < t_value) t = left_child;

  if(left_child + 1 < heap_size) {
    auto right_child = left_child + 1;
    auto right_value = this->contents_[right_child]->Score();
    if(right_value < left_value )
      t = right_child;
  }
  if (pos != t)
  {
    HeapSwap(pos, t);
    Down(t);
  }
}

void GlobalQueue::Up(size_t pos)
{
  while ( this->contents_[pos]->Score()  < this->contents_[pos / 2]->Score() )
  {
    HeapSwap(pos, pos / 2);
    pos >>= 1;
    if(pos==0)
      return;
  }
}

void GlobalQueue::UpdateNode(std::unique_ptr<QueueElement> new_element,size_t pos)
{
  this->contents_[pos] = std::move(new_element);
  Up(pos);
  Down(pos);
}

void GlobalQueue::Insert(std::unique_ptr<QueueElement> queue_element)
{
  auto pos = this->contents_.size();
  // line 26 in algo.5
  queue_element->self_dpp_list_->GQ_pos_ = pos;
  this->contents_.push_back(move(queue_element));
  Up(pos);
}

std::unique_ptr<QueueElement> GlobalQueue::ExtractMin() {
  HeapSwap(0,contents_.size()-1);
  auto top = std::move(contents_[contents_.size()-1]);
  top->self_dpp_list_->GQ_pos_ = static_cast<size_t>(-1);
  contents_.pop_back();
  if(!this->Empty())
    Down(0);
  return top;
}

std::unique_ptr<QueueElement> GlobalQueue::Pop(size_t pos)
{
  HeapSwap(pos,contents_.size()-1);
  auto popped = std::move(contents_[contents_.size()-1]);
  popped->self_dpp_list_->GQ_pos_ = static_cast<size_t>(-1);
  contents_.pop_back();
  Down(pos);
  return popped;
}

