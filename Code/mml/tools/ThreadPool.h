#if !defined MML_THREAD_POOL_H
#define MML_THREAD_POOL_H

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>

namespace MML
{
	class ThreadPool {
	private:
		std::vector<std::thread> workers;
		std::queue<std::function<void()>> tasks;
		std::mutex queue_mutex;
		std::condition_variable condition;
		bool stop;

	public:
		ThreadPool(size_t numThreads) : stop(false)
		{
			for (size_t i = 0; i < numThreads; ++i)
				workers.emplace_back([this]
					{
						while (true) {
							std::function<void()> task;
							{
								std::unique_lock<std::mutex> lock(this->queue_mutex);

								this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });

								if (this->stop && this->tasks.empty())
									return;

								task = std::move(this->tasks.front());
								this->tasks.pop();
							}
							task();
						}
					});
		}
		~ThreadPool()
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			stop = true;
			condition.notify_all();
			
			for (std::thread& worker : workers)
				worker.join();
		}

		// Delete copy constructor and assignment operator
		ThreadPool(const ThreadPool&) = delete;
		ThreadPool& operator=(const ThreadPool&) = delete;

		// Allow move semantics
		ThreadPool(ThreadPool&&) = default;
		ThreadPool& operator=(ThreadPool&&) = default;

		void enqueue(std::function<void()> f) 
		{
			std::unique_lock<std::mutex> lock(queue_mutex);
			tasks.push(std::move(f));

			condition.notify_one();
		}
		void wait_for_tasks() {
			std::unique_lock<std::mutex> lock(queue_mutex);
			condition.wait(lock, [this] { return tasks.empty(); });
		}
		bool has_tasks() {
			std::unique_lock<std::mutex> lock(queue_mutex);
			return !tasks.empty();
		}

	};
} // namespace MML

#endif // MML_THREAD_POOL_H