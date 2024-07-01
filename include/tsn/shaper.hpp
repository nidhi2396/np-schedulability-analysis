#ifndef SHAPER_HPP
#define SHAPER_HPP

#include <ostream>
#include <vector>
#include <deque>
#include <algorithm> // for find
#include <functional> // for hash
#include <exception>
#include <cmath>

#include "time.hpp"
#include "config.h"


namespace NP {

	template<class Time> class Time_Aware_Shaper {

	public:
		typedef std::vector<Time_Aware_Shaper<Time>> TAS_set;
		typedef std::deque<Interval<Time>> Intervals;
		typedef Time Priority; // Make it a time value to support EDF
		typedef Time Period;

	private:
		Priority priority;
		Period period;
		bool TAS;
		bool CBS;
		bool isVar;
		Time guard_band;
		Intervals tas_close_queue;
		mutable unsigned int last_search_result;

	public:

		Time_Aware_Shaper(Priority prio,
			Period per,
			bool tas,
			bool cbs,
			bool isvar,
			Time gb,
			Intervals tascq)
		: priority(prio), period(per), TAS(tas), CBS(cbs), isVar(isvar), guard_band(gb), tas_close_queue(tascq)
		{
			last_search_result=0;
		}

		Time get_period() const
		{
			return period;
		}

		Time get_priority() const
		{
			return priority;
		}

		Time get_guardband() const
		{
			return guard_band;
		}

		bool is_variable() const
		{
			return isVar;
		}

		Time next_open(Time check, Time gband) const
		{
			Time result;
			Intervals gate_close = get_gates_close(check,check,gband);
			result = check;
			for(Interval<Time> gc: gate_close)
			{
				if(check < gc.from()) {
					result = check;
					break;
				} else if(check >= gc.from() && check <= gc.upto()) {
					result = gc.upto() + Time_model::constants<Time>::epsilon();
					break;
				}
			}
			return result;
		}

		Time next_close(Time check) const
		{
			Time result;
			result = Time_model::constants<Time>::infinity();
			unsigned int idx = search_close_queue(check);
			for (auto it = tas_close_queue.begin()+idx; it!=tas_close_queue.end(); it++)
			{
				Interval<Time> gc = *it;
				if (gc.from() == gc.until())
					continue;
				if (check <= gc.from()) {
					result = gc.from();
					break;
				}
			}
			return result;
		}

		Intervals get_gates_open(Time start, Time end, Time gband) const
		{
			Intervals GO;
			if (start > end) {
				return GO;
			}
			Intervals GC = get_gates_close(start, end, gband);
			DM("GC:");
			for(auto st:GC)
				DM(st<<",");
			DM(std::endl);
			if(GC.size() == 0)
			{
				GO.emplace_back(Interval<Time>{start,end});
				return GO;
			}
			Time start_period = floor(start/period);
			Time end_period = ceil(end/period) + 1;
			Time next_gate_open;
			Time first = 0;

			bool skip = false;
			if(start_period*period < GC[0].from())
				next_gate_open = start_period*period;
			else
			{
				next_gate_open = GC[0].upto()+1;
				skip = true;
			}

			for(Interval<Time> gc : GC)
			{
				if(skip == true && first == 0)
				{
					first += 1;
					continue;
				}
				first += 1;
				GO.emplace_back(Interval<Time>{next_gate_open, gc.from() - Time_model::constants<Time>::epsilon()});
				next_gate_open = gc.upto() + Time_model::constants<Time>::epsilon();
			}

			if(next_gate_open < end_period*period)
				GO.emplace_back(Interval<Time>{next_gate_open, end_period*period});

			return GO;
		}

		// The tas_close_queue contains a list of gates. Since the list is
		// likely to cover the complete (hyper) period of the analysis,
		// it can be quite long (10k - 100k).
		// Searching in that list should be very efficient.
		// Binary reseach would still take 15 steps).
		// A hash function might help, if the overhead is not too large.
		//
		// The usage pattern of the search function is that the argument
		// is on average slowly increasing and the return value is typically
		// within 3 of the previous return value. Using the previous value
		// as a starting point might be useful.

		// The search_close_queue(A) function returns the index of the first
		// interval that contain A or is after A.
	  
		unsigned int search_close_queue(Time start) const
		{
			unsigned int idx = last_search_result;
			if (tas_close_queue[idx].upto() < start) {
				// start after the current interval, search forward.
				while (idx<tas_close_queue.size() && tas_close_queue[idx].upto() < start) idx++; 
			} else if (tas_close_queue[idx].from() > start) {
				// start before the current interval, search backward if possible.
				while (idx>0 && tas_close_queue[idx-1].upto() >= start) idx--;
			}
			return idx;
		}

		Intervals get_gates_close(Time start, Time end, Time gband) const
		{
			Intervals mGC;
			if (start > end)
				return mGC;

			Intervals rGC, GC;
			Time start_period = floor(start/period);
			Time end_period = ceil(end/period) + 1;
			Time last_lb=0;
			Time last_ub=0;
			unsigned int idx = search_close_queue(start-(start_period*period));
			
			for(Time current_period = start_period; current_period <= end_period; current_period += 1)
			{
				Time cp = current_period*period;
				for(auto it = tas_close_queue.begin()+idx; it!=tas_close_queue.end(); it++)
				{
					const Interval<Time> gc = *it;
					if (gc.from() == gc.upto())
						continue;

					Time ub = cp + gc.upto() - Time_model::constants<Time>::epsilon();
					if (ub < start)
						continue;

					Time lb = cp + gc.from() - gband;
					if (lb > end)
						break;
					if (lb<ub) {
						if (lb<=last_ub) {
							// Extend.
							last_ub=ub;
						} else {
							// A gab with the previous
							if (last_lb<last_ub)
								mGC.emplace_back(Interval<Time>{last_lb, last_ub});
							last_lb = lb;
							last_ub = ub;
						}
					}
					if (ub > end)
						break;
				}
				idx=0;
			}
			if (last_lb<last_ub) {
				mGC.emplace_back(Interval<Time>{last_lb,last_ub});
			}

			if(mGC.size()==0) {
				return GC;
			}
			// Merge should not be necessary.
			//rGC = merge_ints(mGC);
			// Empty shouldn't occur.
			//GC = remove_empty_ints(rGC);

			return mGC;
		}

		Intervals remove_empty_ints(Intervals remove_lists) const
		{
			Intervals result;
			for(int i=0; i<remove_lists.size();i++)
			{
				if(remove_lists[i].from() == remove_lists[i].upto())
					continue;
				else
					result.emplace_back(remove_lists[i]);
			}
			return result;
		}

		Intervals merge_ints(Intervals merge_lists) const
		{
			Intervals merged;
			Interval<Time> to_add = merge_lists[0];
			for(int i=1;i<merge_lists.size();i++)
			{
				if(to_add.upto() >= merge_lists[i].from())
				{
					to_add = Interval<Time>{to_add.from(),merge_lists[i].upto()};
				}
				else
				{
					merged.emplace_back(to_add);
					to_add = merge_lists[i];
				}
			}
			merged.emplace_back(to_add);
			return merged;
		}
	};

	class InvalidTASParameter : public std::exception
	{
		public:

		InvalidTASParameter(const JobID& bad_id)
		: ref(bad_id)
		{}

		const JobID ref;

		virtual const char* what() const noexcept override
		{
			return "invalid TAS parameter";
		}

	};

	template<class Time>
	void validate_tas_refs(const typename Time_Aware_Shaper<Time>::TAS_set tasQueues,
	                         const typename Job<Time>::Job_set jobs)
	{
		if (tasQueues.size()>0)
		{
			for (auto job : jobs) {
				bool priority_check = false;
				for (auto& tas : tasQueues) {
					if (job.get_priority() == tas.get_priority()) {
						priority_check = true;
						break;
					}
				}
				if (priority_check == false)
					throw InvalidTASParameter(job.get_id());
			}
		}
	}

}

#endif
