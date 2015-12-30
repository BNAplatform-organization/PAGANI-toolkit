#ifndef MAINWINDOW_H
#define MAINWINDOW_H
/*
#include <QMainWindow>
#include <iostream>
#include <QCheckBox>*/

#include "smainwindow.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT   //用来申明信号和槽，因为这是c++所没有的，所以要特别标示，这是一条宏语句

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
     sMainWindow w1;

private slots:
    void on_pushButtonSave_clicked();
    void on_pushButtonLoad_clicked();
 //   void on_toolButtonCUCorMat_Dir_for_BOLD_clicked();
  //  void on_toolButtonLp_input_dir_clicked();
  //  void on_toolButtonCUBFS_Lp_input_dir_clicked();
  //  void on_toolButtonBFS_MulCPU_input_dir_clicked();
//    void on_toolButtonCp_input_dir_clicked();
//    void on_toolButtonDegree_input_dir_clicked();
//    void on_toolButtonCUBC_input_dir_clicked();
//    void on_toolButtonCUEC_input_dir_clicked();
//    void on_toolButtonConvertNII_input_file_clicked();
 //   void on_toolButtonConvertNII_mask_file_clicked();
//    void on_toolButtonL_Modularity_dir_for_csr_clicked();
    void on_toolButtonSaveDir_clicked();
  //  void on_switchButton_clicked();
//    void on_toolButtonPC_CPU_clicked();
    void on_checkCUCorMat_clicked(bool checked);
  // void on_checkLp_clicked(bool checked);
 //   void on_checkCUBFS_Lp_clicked(bool checked);
  //  void on_checkBFS_MulCPU_clicked(bool checked);
    void on_checkCP_clicked(bool checked);
//    void on_checkDegree_clicked(bool checked);
//    void on_checkCUBC_clicked(bool checked);
//    void on_checkCUEC_clicked(bool checked);
    void on_checkL_Modularity_clicked(bool checked);
  //  void on_checkPC_CPU_clicked(bool checked);
   // void on_checkConvertNII_clicked(bool checked);
    //void on_clearscreen();
    void clearscreen();
    void on_toolButton_Working_Directory_clicked();
    void on_toolButton_Mask_File_clicked();

//    void on_toolButtonLp_input_dir_NodalMetrics_clicked();

//    void on_toolButtonCp_input_dir_NodalMetrics_clicked();

//    void on_checkLp_NodalMetrics_clicked(bool checked);

//    void on_checkCP_NodalMetrics_clicked(bool checked);

//    void on_radioButton_CUCorMat_to_average_clicked(bool checked);

    void on_ordinary_clicked(bool checked);

    void on_fisher_clicked(bool checked);

    void on_radioButtonCUCorMat_to_save_cormatrix_clicked(bool checked);

    void on_CheckCUCorMat_to_average_clicked(bool checked);

    void on_switchButton_currentIndexChanged(int index);

   // void on_checkSmallWordProperty_clicked(bool checked);

signals:
    void mySignalMtoS();
    void mySignalMclr();
    void mySignalMgrayBoxCUCorMat(bool);
  //  void mySignalMgrayBoxLp(bool);
  //  void mySignalMgrayBoxCUBFS_Lp(bool);
  //  void mySignalMgrayBoxBFS_MulCPU(bool);
    void mySignalMgrayBoxCP(bool);
    void mySignalMgrayBoxDegree(bool);
    void mySignalMgrayBoxCUBC(bool);
    void mySignalMgrayBoxCUEC(bool);
    void mySignalMgrayBoxL_Modularity(bool);
    void mySignalMgrayBoxPC_CPU(bool);
    void mySignalMgrayBoxConvertNII(bool);
    void mySignalMgrayBoxSmallWorld(bool);
 //   void mySignalMgrayBoxLp_NodalMetrics(bool);
 //   void mySignalMgrayBoxCP_NodalMetrics(bool);

};

#endif // MAINWINDOW_H
