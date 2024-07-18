#pragma once
#include "imageProcessor/logging.h"
#include "containers/ArrayView.h"
#include "containers/ArraySlice.h"
#include "containers/device_array.hpp"

using namespace CUDA;

template<typename T>
class DeviceBufferArray {
public:
	/**
	 * @brief 无参构造函数
	 */
	explicit DeviceBufferArray() : m_buffer(nullptr, 0), m_array(nullptr, 0) {}

	/**
	 * @brief 有参构造函数，指定缓存大小
	 */
	explicit DeviceBufferArray(size_t capacity) {
		AllocateBuffer(capacity);
		m_array = DeviceArray<T>(m_buffer.ptr(), 0);
	}

	/**
	 * @brief 析构函数
	 */
	~DeviceBufferArray() = default;

	//访问与获取的方法
	CUDA::DeviceArray<T> Array() const { return m_array; }
	DeviceArrayView<T> ArrayView() const { return DeviceArrayView<T>(m_array.ptr(), m_array.size()); }
	DeviceArrayView<T> ArrayReadOnly() const { return DeviceArrayView<T>(m_array.ptr(), m_array.size()); }
	DeviceArraySlice<T> ArraySlice() { return DeviceArraySlice<T>(m_array.ptr(), m_array.size()); }
	CUDA::DeviceArray<T> Buffer() const { return m_buffer; }
	
	/**
	 * @brief 与传入的数据交换
	 */
	void swap(DeviceBufferArray<float>& other) {
		m_buffer.swap(other.m_buffer);
		m_array.swap(other.m_array);
	}
	
	//Cast to raw pointer
	const T* Ptr() const { return m_buffer.ptr(); }
	T* Ptr() { return m_buffer.ptr(); }
	operator T*() { return m_buffer.ptr(); }
	operator const T*() const { return m_buffer.ptr(); }
	
	//查询大小
	size_t Capacity() const { return m_buffer.size(); }
	size_t BufferSize() const { return m_buffer.size(); }
	size_t ArraySize() const { return m_array.size(); }

	/**
	 * @brief 分配缓存
	 * 
	 */
	void AllocateBuffer(size_t capacity) {
		if(m_buffer.size() > capacity) return;
		m_buffer.create(capacity);
		m_array = DeviceArray<T>(m_buffer.ptr(), 0);
	}
	/**
	 * @brief 释放缓存
	 * 
	 */
	void ReleaseBuffer() {
		if(m_buffer.size() > 0) m_buffer.release();
	}
	
	/**
	 * @brief 修改数组大小
	 */
	bool ResizeArray(size_t size, bool allocate = false) {
		if(size <= m_buffer.size()) {
			m_array = DeviceArray<T>(m_buffer.ptr(), size);
			return true;
		} 
		else if(allocate) {
			const size_t prev_size = m_array.size();

			//需要先拷贝以前的数据
			DeviceArray<T> old_buffer = m_buffer;
			//分配新的缓存 是要求size的1.5倍
			m_buffer.create(static_cast<size_t>(size * 1.5));
			if(prev_size > 0) {
				cudaSafeCall(cudaMemcpy(m_buffer.ptr(), old_buffer.ptr(), sizeof(T) * prev_size, cudaMemcpyDeviceToDevice));
				old_buffer.release();
			}

			//修改数组大小
			m_array = DeviceArray<T>(m_buffer.ptr(), size);
			return true;
		} 
		else {
			return false;
		}
	}

	/**
	 * @brief 修改数组大小，如果不够则抛出异常
	 */
	void ResizeArrayOrException(size_t size) {
		if (size > m_buffer.size()) {
			LOG(FATAL) << "The pre-allocated buffer is not enough";
		}

		//Change the size of array
		m_array = DeviceArray<T>(m_buffer.ptr(), size);
	}
private:
	DeviceArray<T> m_buffer;
	DeviceArray<T> m_array;
};

