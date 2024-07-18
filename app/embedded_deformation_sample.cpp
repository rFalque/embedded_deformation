#include "imageProcessor/imageProcessor.h"
#include "utils/IO_libIGL/readOBJ.h"
#include "utils/IO/readPLY.h"
#include "utils/IO/writePLY.h"
#include "utils/IO/readCSV.h"



#include "libGraphCpp/readGraphOBJ.hpp"
#include "libGraphCpp/polyscopeWrapper.hpp"

#include "embedded_deformation/embedDeform.hpp"
#include "utils/visualization/plotMesh.h"
#include "utils/visualization/plotCloud.h"
#include "embedded_deformation/options.hpp"

#include <yaml-cpp/yaml.h>
#include "polyscope/polyscope.h"
#include "rigidSolver/rigidSolver.h"

void vertex2image(const Eigen::MatrixXd& source, cv::Mat& result, const Intrinsic& intrinsic)
{
    for(int i = 0; i < source.rows(); i++) {
        Eigen::Vector3d source_point = source.row(i);
        int x = (int) (source_point(0)/(source_point(2) + 1e-10) * intrinsic.focal_x) + intrinsic.principal_x;
        int y = (int) (source_point(1)/(source_point(2) + 1e-10) * intrinsic.focal_y) + intrinsic.principal_y;
        if(x >= 0 && x < result.cols && y >= 0 && y < result.rows) {
            result.at<cv::Vec3d>(y, x) = cv::Vec3d(source_point(0), source_point(1), source_point(2));
        }
    }
}

void diff(const cv::Mat& src_image, const cv::Mat& tag_image, cv::Mat& diff_image)
{
    for(int i = 0; i < src_image.rows; i++) {
        for (int j = 0; j < src_image.cols; j++) {
            cv::Vec3d src = src_image.at<cv::Vec3d>(i,j);
            cv::Vec3d tag = tag_image.at<cv::Vec3d>(i,j);
            double diff = sqrt(pow(src[0] - tag[0], 2) + pow(src[1] - tag[1], 2) + pow(src[2] - tag[2], 2));
            diff_image.at<double>(i,j) = diff;
        }
    }
}

void visualizion(cv::Mat& diff_image1, cv::Mat& diff_image2)
{
    double minVal1, maxVal1, minVal2, maxVal2;
    cv::Point minLoc1, maxLoc1, minLoc2, maxLoc2;
    cv::minMaxLoc(diff_image1, &minVal1, &maxVal1, &minLoc1, &maxLoc1);
    cv::minMaxLoc(diff_image2, &minVal2, &maxVal2, &minLoc2, &maxLoc2);
    double maxVal = max(maxVal1, maxVal2);
    double minVal = min(minVal1, minVal2);
    // 对两个深度图进行归一化 范围是0~255
    cv::Mat diff_image1_normalized = cv::Mat::zeros(diff_image1.rows, diff_image1.cols, CV_8UC3);
    cv::Mat diff_image2_normalized = cv::Mat::zeros(diff_image2.rows, diff_image2.cols, CV_8UC3);
    for( int i = 0; i < diff_image1.rows; i++ )
    {
        for( int j = 0; j < diff_image1.cols; j++ )
        {
            uchar r_value1 = 255*(diff_image1.at<double>(i,j) - minVal)/(maxVal - minVal);
            uchar b_value1 = 255*(maxVal - diff_image1.at<double>(i,j))/(maxVal - minVal);
            diff_image1_normalized.at<cv::Vec3b>(i,j) = static_cast<cv::Vec3b>(r_value1, 
                                                                                0,
                                                                                b_value1);
            diff_image2_normalized.at<cv::Vec3b>(i,j) = static_cast<cv::Vec3b>(255*(diff_image2.at<double>(i,j) - minVal)/(maxVal - minVal), 
                                                                                0,
                                                                                255*(maxVal - diff_image2.at<double>(i,j))/(maxVal - minVal));
        }
    }
    cv::namedWindow("prev&curr", cv::WINDOW_AUTOSIZE);
    cv::namedWindow("deformed&curr", cv::WINDOW_AUTOSIZE);
    cv::imshow("prev&curr", diff_image1_normalized);
    cv::imshow("deformed&curr", diff_image2_normalized);
    cv::waitKey(0);
}

void evaluate(const Eigen::MatrixXd& previous, 
            const Eigen::MatrixXd& previous_deformed, 
            const Eigen::MatrixXd& current, 
            const Intrinsic& intrinsic,
            const unsigned width,
            const unsigned height)
{
    cv::Mat previous_image = cv::Mat::zeros(height, width, CV_64FC3);
    cv::Mat previous_deformed_image = cv::Mat::zeros(height, width, CV_64FC3);
    cv::Mat current_image = cv::Mat::zeros(height, width, CV_64FC3);

    vertex2image(previous, previous_image, intrinsic);
    vertex2image(previous_deformed, previous_deformed_image, intrinsic);
    vertex2image(current, current_image, intrinsic);

    cv::Mat diff_image1 = cv::Mat::zeros(height, width, CV_64FC1);
    cv::Mat diff_image2 = cv::Mat::zeros(height, width, CV_64FC1);

    diff(previous_image, current_image, diff_image1);
    diff(previous_deformed_image, current_image, diff_image2);

    visualizion(diff_image1, diff_image2);
}



int main(int argc, char* argv[])
{
   
    options opts;
    opts.loadYAML("../config.yaml");
    std::cout << "Progress: yaml loaded\n";

    polyscope::init();
    std::cout << "Progress: plyscope initialized\n";

    // 读取图像的类
    std::shared_ptr<FetchInterface> image_fetcher = std::make_shared<FetchInterface>(opts.image_path);

    Intrinsic raw_Intrinsic(opts.focal_x, opts.focal_y, opts.principal_x, opts.principal_y);

    auto image_processor = std::make_shared<imageProcessor>(
                                raw_Intrinsic,                                      
                                opts.width, opts.height, 
                                opts.start_frame,
                                image_fetcher);
    
    // 并行化icp刚性求解
    // To do: 加入dense color icp
    std::shared_ptr<RigidSolver> rigid_solver = std::make_shared<RigidSolver>(image_processor->clip_width(), 
                                                                            image_processor->clip_height(), 
                                                                            image_processor->clip_intrinsic());
    
    // 图像处理 读取图像，计算法向量，计算顶点图，双边滤波，读取有效顶点与法向量
    image_processor->test();
    Eigen::Matrix<float, 3, 4> init_guess = rigid_solver->test(image_processor->getVertexNormalMaps(),image_processor->getVertexNormalMapsPrev());
    std::vector<Surfel> surfels =  image_processor->getSurfelData();
    std::vector<Surfel> surfels_prev = image_processor->getSurfelDataPrev();
    std::vector<correspondence> correspondences = rigid_solver->getCorrespondence();
    
    /* 
     * V: vertex of the surface
     * F: faces of the surface
     * N: nodes of the deformation graph
     * E: edges of the deformation graph
     */
    
    Eigen::MatrixXd V, N;                            // 当前帧顶点图和节点
    Eigen::MatrixXd V_prev, N_prev;                  // 上一帧顶点图和节点
    Eigen::MatrixXi F, E;
    Eigen::MatrixXd normals, normals_prev;           // 当前帧和上一帧的法向量
    EmbeddedDeformation* non_rigid_deformation;
    readPLY("/home/maihn/dev/embeddedDeformation/datasets/breathe_pcd/ply_1.ply", V, F);
    auto* psCloud = polyscope::registerPointCloud("Point Cloud", V);
    psCloud->setPointColor(glm::vec3(0.0, 1.0, 0.0)); // RGB颜色，这里是绿色
    polyscope::getPointCloud("Point Cloud")->setPointRadius(0.00125, true);
    polyscope::show();
    V.resize(surfels.size(), 3);
    for(int i=0; i<surfels.size(); i++)
    {   
        V(i,0) = surfels[i].vertex.x;
        V(i,1) = surfels[i].vertex.y;
        V(i,2) = surfels[i].vertex.z;
    }
    
    
    V_prev.resize(surfels_prev.size(), 3);
    normals_prev.resize(surfels_prev.size(), 3);
    for(int i=0; i<surfels_prev.size(); i++)
    {   
        V_prev(i,0) = surfels_prev[i].vertex.x - 0.5;
        V_prev(i,1) = surfels_prev[i].vertex.y;
        V_prev(i,2) = surfels_prev[i].vertex.z;
        normals_prev(i,0) = surfels_prev[i].normal.x;
        normals_prev(i,1) = surfels_prev[i].normal.y;
        normals_prev(i,2) = surfels_prev[i].normal.z;
    }

    auto* psCloud4 = polyscope::registerPointCloud("Point Cloud 2", V_prev);
    psCloud4->setPointColor(glm::vec3(1.0, 0.0, 0.0)); // RGB颜色，这里是绿色
    polyscope::getPointCloud("Point Cloud 2")->setPointRadius(0.00125, true);
    std::cout<<"using knn distance"<<std::endl;
    
    // 输入的顶点与法向量是所有有效的顶点以及对应的法向量，是上一帧的顶点与法向量
    non_rigid_deformation = new EmbeddedDeformation(V_prev, normals_prev, opts);

    // 展示节点
    // if (opts.visualization)
    //     non_rigid_deformation->show_deformation_graph();
    
    // read correspondences
    // 具有对应关系的点对，相比有效的顶点与法向量，数量会少一些
    Eigen::MatrixXd target_points(correspondences.size(), 3);
    Eigen::MatrixXd source_points(correspondences.size(), 3);
    Eigen::MatrixXd target_normals(correspondences.size(), 3);
    for(int i=0; i<correspondences.size(); i++)
    {
        target_points(i,0) = correspondences[i].tag_vertex.x;
        target_points(i,1) = correspondences[i].tag_vertex.y;
        target_points(i,2) = correspondences[i].tag_vertex.z;
        source_points(i,0) = correspondences[i].src_vertex.x;
        source_points(i,1) = correspondences[i].src_vertex.y;
        source_points(i,2) = correspondences[i].src_vertex.z;
        target_normals(i,0) = correspondences[i].tag_normal.x;
        target_normals(i,1) = correspondences[i].tag_normal.y;
        target_normals(i,2) = correspondences[i].tag_normal.z;
    }

    std::cout << "progress : start deformation ..." << std::endl;
    Eigen::MatrixXd V_deformed, normal_deformed;                       // 变形后的顶点和法向量 
    // 输入的是具有有效对应关系的点对
    std::cout << "progress: non rigid tracking: 1 ..."<<std::endl;
    non_rigid_deformation->deform(source_points, target_points, target_normals, V_deformed, normal_deformed, opts);

    for( int i =1; i<10; i++)
    {
        std::cout<<"progress: non rigid tracking "<< i+1 <<" ..."<<std::endl;
        // 一次非刚性求解后寻找新的匹配点
        DeviceBufferArray<correspondence> new_correspondences;
        int num_pixels = image_processor->clip_width() * image_processor->clip_height();
        new_correspondences.AllocateBuffer(num_pixels);
        rigid_solver->findCorrespondences(V_deformed, normal_deformed, new_correspondences);
        DeviceArrayView<correspondence> new_correspondence_view = new_correspondences.ArrayView();
        std::vector<correspondence> h_correspondences;
        // 将新的匹配点下载到主机
        new_correspondence_view.Download(h_correspondences);
        Eigen::MatrixXd new_target_points(h_correspondences.size(), 3);
        Eigen::MatrixXd new_source_points(h_correspondences.size(), 3);
        Eigen::MatrixXd new_target_normals(h_correspondences.size(), 3);
        for(int i=0; i<h_correspondences.size(); i++)
        {
            new_target_points(i,0) = h_correspondences[i].tag_vertex.x;
            new_target_points(i,1) = h_correspondences[i].tag_vertex.y;
            new_target_points(i,2) = h_correspondences[i].tag_vertex.z;
            new_source_points(i,0) = h_correspondences[i].src_vertex.x;
            new_source_points(i,1) = h_correspondences[i].src_vertex.y;
            new_source_points(i,2) = h_correspondences[i].src_vertex.z;
            new_target_normals(i,0) = h_correspondences[i].tag_normal.x;
            new_target_normals(i,1) = h_correspondences[i].tag_normal.y;
            new_target_normals(i,2) = h_correspondences[i].tag_normal.z;
        }
        V_deformed.resize(0,0);
        normal_deformed.resize(0,0);
        non_rigid_deformation->deform(new_source_points, new_target_points, new_target_normals, V_deformed, normal_deformed, opts);
    }
    
    polyscope::init();
    polyscope::view::upDir = polyscope::view::UpDir::NegYUp;

    // 上一帧 绿色
    auto* psCloud1 = polyscope::registerPointCloud("Point Cloud Previous", V_prev);
    psCloud1->setPointColor(glm::vec3(0.0, 1.0, 0.0)); // RGB颜色，这里是绿色
    // polyscope::getPointCloud("Point Cloud 2")->addVectorQuantity("normals 2", normals_prev);
    polyscope::getPointCloud("Point Cloud Previous")->setPointRadius(0.00125, true);

    // 当前帧 红色
    auto* psCloud2 = polyscope::registerPointCloud("Point Cloud Current", V);
    // polyscope::getPointCloud("Point Cloud 1")->addVectorQuantity("normals 1", normals);
    psCloud2->setPointColor(glm::vec3(1.0, 0.0, 0.0)); // RGB颜色，这里是红色
    polyscope::getPointCloud("Point Cloud Current")->setPointRadius(0.00125, true);

    // 变形后的上一帧 蓝色
    auto* psCloud3 = polyscope::registerPointCloud("Point Cloud Previous Deformed", V_deformed);
    psCloud3->setPointColor(glm::vec3(0.0, 0.0, 1.0)); // RGB颜色，这里是蓝色
    // polyscope::getPointCloud("Point Cloud 2")->addVectorQuantity("normals 2", normals_prev);
    polyscope::getPointCloud("Point Cloud Previous Deformed")->setPointRadius(0.00125, true);
    polyscope::show();


    evaluate(V_prev, V_deformed, V, image_processor->clip_intrinsic(), image_processor->clip_width(), image_processor->clip_height());

    return 0;
}
