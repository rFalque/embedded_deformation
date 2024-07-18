#include "imageProcessor/fetchInterface.h"
#include <opencv2/core.hpp>
#include <boost/filesystem/operations.hpp>
#include <exception>
#include <regex>
#include <iterator>
#include <array>
#include "imageProcessor/logging.h"


const int CV_ANYCOLOR = cv::IMREAD_ANYCOLOR;
const int CV_ANYDEPTH = cv::IMREAD_ANYDEPTH;

FetchInterface::FetchInterface(const std::string data_dir, std::string extension, bool force_no_masks) :
		m_mask_buffer_ix(SIZE_MAX), m_use_masks(false)              // 默认不使用掩码
{
	path data_path(data_dir);
	std::vector<path> sorted_paths;
	if (extension.length() > 0 && extension[0] != '.') {
		extension = "." + extension;
	}
	std::copy(fs::directory_iterator(data_path), fs::directory_iterator(), std::back_inserter(sorted_paths));       // 把data_path 文件下所有目录拷贝到 sorted_path
	std::sort(sorted_paths.begin(), sorted_paths.end());                             								// 默认按照升序排列
	bool unexpected_mask_frame_number = false;
	

	// 把所有的color depth mask 图片的路径放在对应的 vector中
	for (auto& path : sorted_paths) {
		if (FilenameIndicatesDepthImage(path.filename().string(), extension)) {
			int frame_number = GetFrameNumber(path);                                       // 将图片名字中数字字符串转换成整数
			if (frame_number != m_depth_image_paths.size()) {
				throw std::runtime_error("Unexpected depth frame number encountered");
			}
			m_depth_image_paths.push_back(path);
		} else if (FilenameIndicatesRGBImage(path.filename().string(), extension)) {
			int frame_number = GetFrameNumber(path);
			if (frame_number != m_rgb_image_paths.size()) {
				throw std::runtime_error("Unexpected RGB frame number encountered");
			}
			m_rgb_image_paths.push_back(path);
		} else if (!force_no_masks && FilenameIndicatesMaskImage(path.filename().string(), extension)) {
			int frame_number = GetFrameNumber(path);
			if (frame_number != m_mask_image_paths.size()) {
				unexpected_mask_frame_number = true;
			}
			m_mask_image_paths.push_back(path);
		}
	}


	
	// 如果depth和rgb的number不一致
	if (m_depth_image_paths.size() != m_rgb_image_paths.size()) {
		LOG(FATAL) << "Presumed depth image count (" << m_depth_image_paths.size() <<
		           ") doesn't equal presumed rgb image count." << m_rgb_image_paths.size();
	}

	// 是否强制不适用mask 默认是false
	if (!force_no_masks) {
		if (unexpected_mask_frame_number) {
			LOG(WARNING)
					<< "Warning: inconsistent mask frame numbers encountered in the filenames. Proceeding without masks.";
		} else if (!m_mask_image_paths.empty()) {               // 掩码路径不为空
			if (m_depth_image_paths.size() != m_mask_image_paths.size() && !m_mask_image_paths.empty()) {       // mask与depthshu数量不一致
				LOG(WARNING)
						<< "Warning: seems like there were some mask image files, but their number doesn't match the "
						   "number of depth frames. Proceeding without masks.";
			} else {
				m_use_masks = true;
			}
		}
	}
}

bool FetchInterface::HasSubstringFromSet(const std::string& string, const std::string* set, int set_size)
{
	bool found_indicator = false;
	for (int i_target_string = 0; i_target_string < set_size; i_target_string++) {
		if (string.find(set[i_target_string]) != std::string::npos) {                          // 在字符串 string 中查找子字符串
			return true;
		}
	}
	return false;
}

bool FetchInterface::FilenameIndicatesDepthImage(const path& filename, const std::string& valid_extension)
{
	if (filename.extension() != valid_extension) return false;                                                                // 首先检查扩展名
	const std::array<std::string, 3> possible_depth_indicators = {"depth", "DEPTH", "Depth"};
	return HasSubstringFromSet(filename.string(), possible_depth_indicators.data(), possible_depth_indicators.size());
}

bool FetchInterface::FilenameIndicatesRGBImage(const path& filename, const std::string& valid_extension)
{
	if (filename.extension() != valid_extension) return false;
	const std::array<std::string, 5> possible_depth_indicators = {"color", "COLOR", "Color", "rgb", "RGB"};
	return HasSubstringFromSet(filename.string(), possible_depth_indicators.data(), possible_depth_indicators.size());
}


bool FetchInterface::FilenameIndicatesMaskImage(const path& filename,
                                                              const std::string& valid_extension)
{
	if (filename.extension() != valid_extension) return false;
	const std::array<std::string, 3> possible_depth_indicators = {"mask", "Mask", "MASK"};
	return HasSubstringFromSet(filename.string(), possible_depth_indicators.data(), possible_depth_indicators.size());
}

int FetchInterface::GetFrameNumber(const path& filename)        // 将图片名字中数字字符串转换成整数
{
	const std::regex digits_regex("\\d+");
	std::smatch match_result;
	const std::string filename_stem = filename.stem().string();
	if (!std::regex_search(filename_stem, match_result, digits_regex)) {
		throw std::runtime_error("Could not find frame number in filename.");
	};
	return std::stoi(match_result.str(0));
}

/**
 * @brief 根据深度图像路径读取深度图
 * 
 */
void FetchInterface::FetchDepthImage(size_t frame_idx, cv::Mat& depth_img)
{
	path file_path = this->m_depth_image_paths[frame_idx];
	// Read the image
	depth_img = cv::imread(file_path.string(), CV_ANYCOLOR | CV_ANYDEPTH);
	// 如果提供了掩码图，那么将深度图中掩码图为0的像素置为0
	if (this->m_use_masks) {
		mask_mutex.lock();
		if (this->m_mask_buffer_ix != frame_idx) {
			m_mask_image_buffer = cv::imread(this->m_mask_image_paths[frame_idx].string(),
											 CV_ANYCOLOR | CV_ANYDEPTH);
			m_mask_buffer_ix = frame_idx;
		}
		cv::Mat masked;
		// m_mask_image_buffer的非零元素对应的depth_img的元素保留，其余元素置为0
		depth_img.copyTo(masked, m_mask_image_buffer);
		mask_mutex.unlock();
		depth_img = masked;
	}
}


void FetchInterface::FetchRGBImage(size_t frame_idx, cv::Mat& rgb_img)
{
	path file_path = this->m_rgb_image_paths[frame_idx];
	//Read the image
	rgb_img = cv::imread(file_path.string(), CV_ANYCOLOR | CV_ANYDEPTH); 

	if (this->m_use_masks) {
		mask_mutex.lock();
		if (this->m_mask_buffer_ix != frame_idx) {
			m_mask_image_buffer = cv::imread(this->m_mask_image_paths[frame_idx].string(),
			                                 CV_ANYCOLOR | CV_ANYDEPTH); 
			m_mask_buffer_ix = frame_idx;
		}
		cv::Mat masked;
		// m_mask_image_buffer的非零元素对应的depth_img的元素保留，其余元素置为0
		rgb_img.copyTo(masked, m_mask_image_buffer);
		mask_mutex.unlock();
		rgb_img = masked;
	}
}
