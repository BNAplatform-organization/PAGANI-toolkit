#ifndef SMAINWINDOW_H
#define SMAINWINDOW_H

#include <QMainWindow>
#include <iostream>
#include <QCheckBox>

#include <QMainWindow>

extern QString file_name_transimit;

extern enum OperatingSystem {
    os_win32,
    os_linux
} operating_system;


namespace Ui {
class sMainWindow;
}

class sMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit sMainWindow(QWidget *parent = 0);
    ~sMainWindow();
   void setParent(QWidget *parent);


private slots:


    void on_pushButtonSave_clicked();

    void on_pushButtonLoad_clicked();

  //  void on_toolButtonCUCorMat_Dir_for_BOLD_clicked();

  //  void on_toolButtonCUBFW_Lp_input_dir_clicked();

 //   void on_toolButtonPC_CPU_input_dir_clicked();

   //  void on_toolButtonCUCP_input_dir_clicked();

  //  void on_toolButtonCp_input_dir_clicked();

 //   void on_toolButtonDegree_input_dir_clicked();

 //   void on_toolButtonCUEC_input_dir_clicked();

//    void on_toolButtonConvertNII_input_file_clicked();

 //   void on_toolButtonConvertNII_mask_file_clicked();

 //   void on_toolButtonL_Modularity_dir_for_csr_clicked();

    void on_toolButtonSaveDir_clicked();

    void on_switchButton_clicked();

    void on_checkCUCorMat_clicked(bool checked);

    void on_checkCUBFW_Lp_clicked(bool checked);

 //   void on_checkCUCP_clicked(bool checked);

    void on_checkCP_clicked(bool checked);

//    void on_checkDegree_clicked(bool checked);

//    void on_checkCUEC_clicked(bool checked);

    void on_checkL_Modularity_clicked(bool checked);

    void on_checkPC_CPU_clicked(bool checked);

    void on_checkConvertNII_clicked(bool checked);
   void clearscreen();
   void on_toolButton_Working_Directory_clicked();

   void on_toolButton_Mask_File_clicked();

   void on_checkCP_Nodal_Metrics_clicked(bool checked);

   void on_radioButton_CUCorMat_to_average_clicked(bool checked);

   void on_fisher_clicked(bool checked);

   void on_ordinary_clicked(bool checked);

   void on_radioButtonCUCorMat_to_save_cormatrix_clicked(bool checked);

signals:
     void mySignalStoM();
     void mySignalSclr();
     void mySignalSgrayBoxCUCorMat(bool);
     void mySignalSgrayBoxCUBFW_Lp(bool);
     void mySignalSgrayBoxCUEC(bool);
     void mySignalSgrayBoxCP(bool);
     void mySignalSgrayBoxDegree(bool);
     void mySignalSgrayBoxCUCP(bool);
     void mySignalSgrayBoxL_Modularity(bool);
     void mySignalSgrayBoxPC_CPU(bool);
     void mySignalSgrayBoxConvertNII(bool);
     void mySignalSgrayBoxCP_Nodal_Metrics(bool);
    // void mySignalSgrayBoxCUCorMat(B);
private:
    Ui::sMainWindow *ui;
     QWidget* pWid;
};

#endif // SMAINWINDOW_H
