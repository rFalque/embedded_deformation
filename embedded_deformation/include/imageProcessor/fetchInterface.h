#pragma once
#include <opencv2/opencv.hpp>
#include <memory>
#include <boost/filesystem/path.hpp>
#include <vector>
#include <mutex>

namespace fs = boost::filesystem;
using path = boost::filesystem::path;

class FetchInterface
{
public:
    
    //Default contruct and de-construct
    FetchInterface(const std::string data_dir, std::string extension = ".png", bool force_no_masks = false);
    ~FetchInterface() = default;

    //Buffer may be maintained outside fetch object for thread safety
    void FetchDepthImage(size_t frame_idx, cv::Mat& depth_img);                  // 虚函数

    //Should be rgb, in CV_8UC3 format
    void FetchRGBImage(size_t frame_idx, cv::Mat& rgb_img);
private:
    static int GetFrameNumber(const path& filename);
    static bool HasSubstringFromSet(const std::string& string, const std::string* set, int set_size);
    static bool FilenameIndicatesDepthImage(const path& filename, const std::string& valid_extension);
    static bool FilenameIndicatesRGBImage(const path& filename, const std::string& valid_extension);
    static bool FilenameIndicatesMaskImage(const path& filename, const std::string& valid_extension);

    std::vector<path> m_rgb_image_paths;
    std::vector<path> m_depth_image_paths;
    std::vector<path> m_mask_image_paths;

    size_t m_mask_buffer_ix;
    cv::Mat m_mask_image_buffer;
    bool m_use_masks;
    std::mutex mask_mutex;
};

