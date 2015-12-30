#include "mainwindow.h"
#include "ui_mainwindow.h"



#include <QFileDialog>
#include <QMessageBox>
#include <QtGlobal>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using  std::string;
extern QString file_name_transimit;
bool both_metrics=0;
bool average_ordinary=0;
bool average_fisher=0;
string to_save_cormatrix="n";
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
 // this->setWindowFlags(windowFlags()| Qt::WindowMaximizeButtonHint);
  // this->setWindowFlags( windowFlags() | Qt::WindowMaximizeButtonHint | Qt::WindowMinimizeButtonHint | Qt::WindowStaysOnTopHint );
    connect(this,SIGNAL(mySignalMtoS()),&w1,SLOT(on_pushButtonLoad_clicked()));
    connect(&w1,SIGNAL(mySignalStoM()),this,SLOT(on_pushButtonLoad_clicked()));
    //这种connect的方式我也是醉了！！！！！
  //  connect(this,SIGNAL(mySignalMclr()),this,SLOT(on_clearscreen()));
    connect(this,SIGNAL(mySignalMclr()),this,SLOT(clearscreen()));


     connect(this,SIGNAL(mySignalMgrayBoxCUCorMat(bool)),this,SLOT(on_checkCUCorMat_clicked(bool)));
  //   connect(this,SIGNAL(mySignalMgrayBoxLp(bool)),this,SLOT(on_checkLp_clicked(bool)));
    // connect(this,SIGNAL(mySignalMgrayBoxCUBFS_Lp(bool)),this,SLOT(on_checkCUBFS_Lp_clicked(bool)));
    // connect(this,SIGNAL(mySignalMgrayBoxBFS_MulCPU(bool)),this,SLOT(on_checkBFS_MulCPU_clicked(bool)));
     connect(this,SIGNAL(mySignalMgrayBoxCP(bool)),this,SLOT(on_checkCP_clicked(bool)));
   //  connect(this,SIGNAL(mySignalMgrayBoxDegree(bool)),this,SLOT(on_checkDegree_clicked(bool)));
   //  connect(this,SIGNAL(mySignalMgrayBoxCUBC(bool)),this,SLOT(on_checkCUBC_clicked(bool)));
   //  connect(this,SIGNAL(mySignalMgrayBoxCUEC(bool)),this,SLOT(on_checkCUEC_clicked(bool)));
     connect(this,SIGNAL(mySignalMgrayBoxL_Modularity(bool)),this,SLOT(on_checkL_Modularity_clicked(bool)));
   //  connect(this,SIGNAL(mySignalMgrayBoxPC_CPU(bool)),this,SLOT(on_checkPC_CPU_clicked(bool)));
  //   connect(this,SIGNAL(mySignalMgrayBoxConvertNII(bool)),this,SLOT(on_checkConvertNII_clicked(bool)));

   //  connect(this,SIGNAL(mySignalMgrayBoxLp_NodalMetrics(bool)),this,SLOT(on_checkLp_NodalMetrics_clicked(bool)));
   //  connect(this,SIGNAL(mySignalMgrayBoxCP_NodalMetrics(bool)),this,SLOT(on_checkCP_NodalMetrics_clicked(bool)));

     connect(this,SIGNAL(mySignalMgrayBoxSmallWorld(bool)),this,SLOT(on_checkSmallWordProperty_clicked(bool)));

         ui->lineEdit_Working_Directory->setText("sdffsssssssss");
     emit mySignalMclr();




}

void MainWindow::clearscreen(){
    qDebug("clear screen!");
    //清屏函数两个功能：1.清空内容,顺带初始化；2.使界面回到灰框状态
    //1.清空内容；
  //  ui->lineEditCUCorMat_Dir_for_BOLD->setText("");
    ui->mask_threshold->setText("");
    ui->lineEdit_Mask_File->setText("");
    ui->CheckCUCorMat_to_average->setEnabled(false);
    ui->ordinary->setChecked(true);
    ui->ordinary->setEnabled(false);
    ui->fisher->setEnabled(false);
   // ui->CheckCUCorMat_to_average->setChecked(false);


    ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(false);
    ui->comboBoxCUCorMat_threshold_type->setCurrentIndex(0);
    ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setText("");

    ui->lineEdit_Working_Directory->setText("");;
    ui->lineEditLp_num_of_random_networks->setText("");

   // ui->comboBoxCp_Cp_type->setCurrentIndex(0);
  //  ui->toolButtonLp_input_dir->setEnabled(false);

    //ui->lineEditCUBFS_Lp_input_dir->setText("");
    //ui->lineEditCUBFS_Lp_num_of_random_networks->setText("");
/*
    ui->lineEditBFS_MulCPU_input_dir->setText("");
    ui->lineEditBFS_MulCPU_num_of_random_networks->setText("");
*/
  //  ui->lineEditCp_input_dir->setText("");
  //  ui->lineEditCp_num_of_random_networks->setText("");



       ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(0);
     //  ui->lineEditL_Modularity_dir_for_csr->setText("");
       ui->lineEditLp_num_of_random_networks->setText("");

       //2.使界面回到灰框状态 good method
       emit mySignalMgrayBoxCUCorMat(false);
    //   emit mySignalMgrayBoxLp(false);
       //emit mySignalMgrayBoxCUBFS_Lp(false);
     //  emit mySignalMgrayBoxBFS_MulCPU(false);
       emit mySignalMgrayBoxCP(false);
     //  emit mySignalMgrayBoxDegree(false);
    //   emit mySignalMgrayBoxCUBC(false);
     //  emit mySignalMgrayBoxCUEC(false);
       emit mySignalMgrayBoxL_Modularity(false);
    //   emit mySignalMgrayBoxPC_CPU(false);
     //  emit mySignalMgrayBoxConvertNII(false);
      // emit mySignalMgrayBoxLp_NodalMetrics(false);
     //  emit mySignalMgrayBoxCP_NodalMetrics(false);
       emit mySignalMgrayBoxSmallWorld(false);

}
/*
void MainWindow::on_clearscreen()
{
    qDebug("what why???");
    ui->mask_threshold->setText("");

    ui->CheckCUCorMat_to_average->setEnabled(false);
    ui->ordinary->setChecked(true);
    ui->ordinary->setEnabled(false);
    ui->fisher->setEnabled(false);
   // ui->CheckCUCorMat_to_average->setChecked(false);


    ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(false);
    ui->comboBoxCUCorMat_threshold_type->setCurrentIndex(0);
    ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setText("");

    ui->lineEdit_Working_Directory->setText("");;
    ui->lineEditLp_num_of_random_networks->setText("");

   // ui->comboBoxCp_Cp_type->setCurrentIndex(0);
  //  ui->toolButtonLp_input_dir->setEnabled(false);

    //ui->lineEditCUBFS_Lp_input_dir->setText("");
    //ui->lineEditCUBFS_Lp_num_of_random_networks->setText("");
/*
    ui->lineEditBFS_MulCPU_input_dir->setText("");
    ui->lineEditBFS_MulCPU_num_of_random_networks->setText("");
*/
  //  ui->lineEditCp_input_dir->setText("");
  //  ui->lineEditCp_num_of_random_networks->setText("");

   //   ui->lineEditDegree_input_dir->setText("");

   //   ui->lineEditCUBC_input_dir->setText("");

    //   ui->lineEditCUEC_input_dir->setText("");

  /*     ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(0);
     //  ui->lineEditL_Modularity_dir_for_csr->setText("");
       ui->lineEditLp_num_of_random_networks->setText("");*/

   //    ui->lineEditPC_CPU_input_dir->setText("");
   //    ui->comboBoxPC_CPU_type_for_participant_coefficient->setCurrentIndex(0);

   //    ui->lineEditConvertNII_input_file->setText("");
   //    ui->lineEditConvertNII_mask_file->setText("");
       //ui->lineEditConvertNII_mask_threshold->setText("");


    //   ui->lineEditLp_input_dir_NodalMetrics->setText("");;
       //ui->lineEditLp_num_of_random_networks_NodalMetrics->setText("");
    //   ui->toolButtonLp_input_dir_NodalMetrics->setEnabled(false);

       //ui->lineEditCUBFS_Lp_input_dir->setText("");
       //ui->lineEditCUBFS_Lp_num_of_random_networks->setText("");
   /*
       ui->lineEditBFS_MulCPU_input_dir->setText("");
       ui->lineEditBFS_MulCPU_num_of_random_networks->setText("");
   */
     //  ui->lineEditCp_input_dir_NodalMetrics->setText("");
    //   ui->lineEditCp_num_of_random_networks_NodalMetrics->setText("");
   /*  qDebug("what why???");
}*/
MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButtonSave_clicked()
{
    if (operating_system != os_win32 && operating_system != os_linux) {
        std::cout << "Not win32 or linux!" << std::endl;

    }
//temporary sequence:lp,cp,smallworld,modularity,pc,degree,bc,ec,outputNII
    std::stringstream script;

  if(ui->switchButton->currentIndex()==0)
  {
     script <<  "echo unweightednetworks" << std::endl;
     string unweighted=ui->lineEdit_Working_Directory->text().toStdString().append("/").append("unweighted");
     if (ui->checkCUCorMat->isChecked()) {
        //组合之外的统统选n
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->mask_threshold->text().isEmpty()
                //|| ui->lineEditCUCorMat_to_average->text().isEmpty()
                //|| ui->lineEditCUCorMat_to_save_cormatrix->text().isEmpty()
                //|| ui->lineEditCUCorMat_threshold_type->text().isEmpty()
                || ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string  type=ui->comboBoxCUCorMat_threshold_type->currentText().toStdString()=="correlation"?"r":"s";
        if(both_metrics==0)
            script << (operating_system == os_win32 ? ".\\exefiles\\CUCorMat.exe " : "./exefiles/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              ' ' <<
              "n" <<
              ' ' <<
              to_save_cormatrix <<
              ' ' <<
              type<<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else if(both_metrics==1&&average_fisher==0){
            script << (operating_system == os_win32 ? ".\\exefiles\\CUCorMat.exe " : "./exefiles/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              ' ' <<
              "bn" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              to_save_cormatrix <<
                  //ui->lineEditCUCorMat_to_save_cormatrix->text().toStdString() <<
              ' ' <<
              type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;}
        //&&average_ordinary==0
        else if(both_metrics==1&&average_fisher==1)
            script << (operating_system == os_win32 ? ".\\exefiles\\CUCorMat.exe " : "./exefiles/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              ' ' <<
              "bf" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              to_save_cormatrix <<
                  //ui->lineEditCUCorMat_to_save_cormatrix->text().toStdString() <<
              ' ' <<
             type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
    }

    if ((ui->checkLp->isChecked())&&(!ui->checkLp_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
          //    ui->lineEdit_Working_Directory->text().toStdString() <<
            unweighted<<
                  ' ' <<
              0 <<' ' <<"g"<<
              std::endl;
    }
    if ((!ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
          //    ui->lineEdit_Working_Directory->text().toStdString() <<
            unweighted<<
                  ' ' <<
              0 <<' ' <<"n"<<
              std::endl;

    }
    if ((!ui->checkLp->isChecked())&&(!ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"l"<<
              std::endl;
    }
    if ((ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
          //    ui->lineEdit_Working_Directory->text().toStdString() <<
            unweighted<<
                  ' ' <<
              0 <<' ' <<"gn"<<
              std::endl;

    }
    if ((ui->checkLp->isChecked())&&(!ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"gl"<<
              std::endl;

    }
    if ((!ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"nl"<<
              std::endl;
    }
    if ((ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"gnl"<<
              std::endl;
    }
        /*
        if (ui->checkLp->isChecked()&&(!ui->checkLp_NodalMetrics->isChecked())) {



    }
    if (ui->checkLp_NodalMetrics->isChecked()&&(!ui->checkLp_NodalMetrics->isChecked())) {

        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"n"<<
              std::endl;
    }
    if (ui->checkLp->isChecked()&&ui->checkLp_NodalMetrics->isChecked()) {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
               unweighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"b"<<
               std::endl;
     }*/
    if (ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
             0 <<' '<<"g"<<
              std::endl;
    }
    if (!ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
             0 <<' '<<"n"<<
              std::endl;
    }
    if (!ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"k"<<
              std::endl;
    }
    if (ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
             0 <<' '<<"gn"<<
              std::endl;
    }
    if (ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"gk"<<
              std::endl;
    }
    if (!ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"nk"<<
              std::endl;
    }
    if (ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
    {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"gnk"<<
              std::endl;
    }

    /*
    if (ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
              unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"g"<<
              std::endl;
    }
    if (ui->checkCP_NodalMetrics->isChecked()&&(!ui->checkCP->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
             unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"n"<<
              std::endl;
    }
    if (ui->checkCP->isChecked()&&ui->checkCP_NodalMetrics->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Cp.exe " : "./exefiles/Cp ") <<
             unweighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"b"<<
              std::endl;
    }*/
    if(ui->checkSmallWordProperty->isChecked()){
        if (ui->lineEdit_Working_Directory->text().isEmpty()||
                ui->lineEditLp_num_of_random_networks->text().isEmpty())
        {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\SmallWorldProperty.exe " : "./exefiles/SmallWorldProperty ") <<
                   unweighted <<
                    ' ' <<
                        std::endl;
    }
    if (ui->checkL_Modularity->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
              //  || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }                                           //这句话是看是不是win32 是就是一种格式 不是就是另一种格式  \\是c++中的转义字符  代表win中的反斜杠；系统路径，编程路径 \\或/ 理解：因为\会被编译错

        string s=ui->comboBoxLouvain_Modularity_modularity_type->currentText().toStdString();
        if(s=="Louvain")
           script << (operating_system == os_win32 ? ".\\exefiles\\Louvain_Modularity.exe " : "./exefiles/Louvain_Modularity ") <<
                 unweighted <<
                 ' ' <<
                 ui->lineEditLp_num_of_random_networks->text().toStdString() <<
                 std::endl;
        else if (s=="Newman_GPU")
            script << (operating_system == os_win32 ? ".\\exefiles\\CUmodularity.exe " : "./exefiles/CUmodularity ") <<
                        unweighted <<
                        ' ' <<
                        ui->lineEditLp_num_of_random_networks->text().toStdString() <<
                        std::endl;
         else
            script << (operating_system == os_win32 ? ".\\exefiles\\Newman_Modularity_CPU.exe " : "./exefiles/Newman_Modularity_CPU ") <<
                        unweighted<<
                         ' ' <<
                         ui->lineEditLp_num_of_random_networks->text().toStdString() <<
                         std::endl;

    }
    if (ui->checkPC_CPU->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
            script << (operating_system == os_win32 ? ".\\exefiles\\PC_CPU.exe " : "./exefiles/PC_CPU ") <<
                        unweighted <<
                        ' ' <<
                            std::endl;
    }
    if (ui->checkDegree->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\Degree.exe " : "./exefiles/Degree ") <<
              unweighted <<
              std::endl;
    }
    if (ui->checkCUBC->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUBC.exe " : "./exefiles/CUBC ") <<
              unweighted <<
              std::endl;
    }
    if (ui->checkCUEC->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\CUEC.exe " : "./exefiles/CUBC ") <<
              unweighted <<
              std::endl;
    }
    if (ui->checkConvertNII->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEdit_Mask_File->text().isEmpty()
                || ui->mask_threshold->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles\\ConvertNII.exe " : "./exefiles/ConvertNII ") <<
              unweighted <<
              ' ' <<
              ui->lineEdit_Mask_File->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              std::endl;
    }
   }
  if(ui->switchButton->currentIndex()==1)
  {
    script <<  "echo weightednetworks" << std::endl;
    string weighted=ui->lineEdit_Working_Directory->text().toStdString().append("/").append("weighted");
     if (ui->checkCUCorMat->isChecked()) {
        //组合之外的统统选n
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->mask_threshold->text().isEmpty()
                //|| ui->lineEditCUCorMat_to_average->text().isEmpty()
                //|| ui->lineEditCUCorMat_to_save_cormatrix->text().isEmpty()
                //|| ui->lineEditCUCorMat_threshold_type->text().isEmpty()
                || ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string  type=ui->comboBoxCUCorMat_threshold_type->currentText().toStdString()=="correlation"?"r":"s";
        if(both_metrics==0)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              ' ' <<
              "n" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              to_save_cormatrix <<
                  //ui->lineEditCUCorMat_to_save_cormatrix->text().toStdString() <<
              ' ' <<
              type<<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else if(both_metrics==1&&average_fisher==0)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              ' ' <<
              "bn" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              to_save_cormatrix <<
                  //ui->lineEditCUCorMat_to_save_cormatrix->text().toStdString() <<
              ' ' <<
              type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else if(both_metrics==1&&average_ordinary==0&&average_fisher==1)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              ' ' <<
              "bf" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              to_save_cormatrix <<
                  //ui->lineEditCUCorMat_to_save_cormatrix->text().toStdString() <<
              ' ' <<
             type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
    }


   /* if (ui->checkLp->isChecked()&&(!ui->checkLp_NodalMetrics->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
              ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
       script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe " : "./exefiles_weighted/CUBFW_Lp ") <<
            weighted<<
             ' ' <<
             ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"g"<<
             std::endl;
    }
    if (ui->checkLp_NodalMetrics->isChecked()&&(!ui->checkLp_NodalMetrics->isChecked())) {

        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe " : "./exefiles_weighted/CUBFW_Lp ") <<
              weighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' '<<"n"<<
              std::endl;
    }
    if (ui->checkLp->isChecked()&&ui->checkLp_NodalMetrics->isChecked()) {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe " : "./exefiles_weighted/CUBFW_Lp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"b"<<
               std::endl;
     }
   */
     if ((ui->checkLp->isChecked())&&(!ui->checkLp_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
           //    ui->lineEdit_Working_Directory->text().toStdString() <<
             weighted<<
                   ' ' <<
               0 <<' ' <<"g"<<
               std::endl;
     }
     if ((!ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
           //    ui->lineEdit_Working_Directory->text().toStdString() <<
             weighted<<
                   ' ' <<
               0 <<' ' <<"n"<<
               std::endl;

     }
     if ((!ui->checkLp->isChecked())&&(!ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"l"<<
               std::endl;
     }
     if ((ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
           //    ui->lineEdit_Working_Directory->text().toStdString() <<
             weighted<<
                   ' ' <<
               0 <<' ' <<"gn"<<
               std::endl;

     }
     if ((ui->checkLp->isChecked())&&(!ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"gl"<<
               std::endl;

     }
     if ((!ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"nl"<<
               std::endl;
     }
     if ((ui->checkLp->isChecked())&&(ui->checkLp_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_lambda->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         script << (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe " : "./exefiles/CUBFW_Lp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<"gnl"<<
               std::endl;
     }


     if (ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_gamma->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }

         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               0 <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"g"<<
               std::endl;
     }
     if (!ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_gamma->isChecked()))
     {

         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }

         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               0 <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"n"<<
               std::endl;
     }
     if (!ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"k"<<
               std::endl;
     }
     if (ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(!ui->checkSmallWordProperty_gamma->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 ) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }

         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               0 <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"gn"<<
               std::endl;
     }
     if (ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"gk"<<
               std::endl;
     }
     if (!ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"nk"<<
               std::endl;
     }
     if (ui->checkCP->isChecked()&&(ui->checkCP_NodalMetrics->isChecked())&&(ui->checkSmallWordProperty_gamma->isChecked()))
     {
         if (ui->lineEdit_Working_Directory->text().isEmpty()
                 || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
             QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         if (ui->lineEditLp_num_of_random_networks->text()=="0")
         {
             QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
             return;
         }
         string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
         script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
               weighted <<
               ' ' <<
               ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"gnk"<<
               std::endl;
     }

     /*if (ui->checkCP->isChecked()&&(!ui->checkCP_NodalMetrics->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
              weighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"g"<<
              std::endl;
    }
    if (ui->checkCP_NodalMetrics->isChecked()&&(!ui->checkCP->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
             weighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"n"<<
              std::endl;
    }
    if (ui->checkCP->isChecked()&&ui->checkCP_NodalMetrics->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditLp_num_of_random_networks->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string s_Cp=ui->comboBoxCp_Cp_type->currentText().toStdString();
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
              weighted <<
              ' ' <<
              ui->lineEditLp_num_of_random_networks->text().toStdString() <<' ' <<(s_Cp == "Onnela" ? "2 " : " 1 ")<<' '<<"b"<<
              std::endl;
    }
   */

    if(ui->checkSmallWordProperty->isChecked()){
        if (ui->lineEdit_Working_Directory->text().isEmpty()||
                ui->lineEditLp_num_of_random_networks->text().isEmpty())
        {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if (ui->lineEditLp_num_of_random_networks->text()=="0")
        {
            QMessageBox::information(this, "Warning", "The value of random networkss(n) can't be zero.", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\SmallWorldProperty.exe " : "./exefiles_weighted/SmallWorldProperty ") <<
                    weighted <<
                    ' ' <<
                        std::endl;
    }
    if (ui->checkL_Modularity->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
               ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }                                           //这句话是看是不是win32 是就是一种格式 不是就是另一种格式  \\是c++中的转义字符  代表win中的反斜杠；系统路径，编程路径 \\或/ 理解：因为\会被编译错

           script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Louvain_Modularity.exe " : "./exefiles_weighted/Louvain_Modularity ") <<
               weighted <<
                 ' ' <<
                 ui->lineEditLp_num_of_random_networks->text().toStdString() <<
                 std::endl;
    }
    if (ui->checkPC_CPU->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
            script << (operating_system == os_win32 ? ".\\exefiles_weighted\\PC_CPU.exe " : "./exefiles_weighted/PC_CPU ") <<
                        weighted <<
                        ' ' <<
                            std::endl;
    }
    if (ui->checkDegree->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Degree.exe " : "./exefiles_weighted/Degree ") <<
             weighted <<
              std::endl;
    }
    if (ui->checkCUBC->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBC.exe " : "./exefiles_weighted/CUBC ") <<
             weighted <<
              std::endl;
    }
    if (ui->checkCUEC->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUEC.exe " : "./exefiles_weighted/CUBC ") <<
              weighted <<
              std::endl;
    }
    if (ui->checkConvertNII->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEdit_Mask_File->text().isEmpty()
                || ui->mask_threshold->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\ConvertNII.exe " : "./exefiles_weighted/ConvertNII ") <<
             weighted <<
              ' ' <<
              ui->lineEdit_Mask_File->text().toStdString() <<
              ' ' <<
              ui->mask_threshold->text().toStdString() <<
              std::endl;
    }
   }

    QString file_name = ui->lineEditSaveDir->text();
    std::cout << script.str() << std::endl;
    QFileInfo file_info(file_name);//是指文件不存在，只有text内容为空文件才会不存在
    if (!file_info.exists())
    {
        file_name = QFileDialog::getSaveFileName(this,
                                                         "Save as..." ,
                                                         (operating_system == os_win32 ? "script.bat" : "script.sh"),
                                                         (operating_system == os_win32 ? "script (*.bat);;Any (*.*)" : "script (*.sh);Any (*.*)"));

    }

    //filename为空，赋值为null
    if (!file_name.isNull()) { //不空1.怎么知道os是哪种输出设备？2.
        ui->lineEditSaveDir->setText(file_name);
        std::ofstream os;
        os.open(file_name.toStdString().c_str());
        os << script.str();//.str表示一份拷贝
        os.close();   //保存完了

        if (operating_system == os_win32) {
            std::string cmd = "start " + file_name.toStdString() + " &";
            system(cmd.c_str());
        } else {
            std::string cmd = "sh " + file_name.toStdString() + " &";
            system(cmd.c_str());
        }    //执行命令
    }
    /*if (ui->lineEditSaveDir->text().isEmpty()) {
        on_toolButtonSaveDir_clicked();
    }
    QString file_name = ui->lineEditSaveDir->text();
    if (!file_name.isNull()) {
        file_name += ((operating_system == win32) ? "/script.bat" : "script.sh");
        std::ofstream os;
        os.open(file_name.toStdString().c_str());
        os << script.str();
        os.close();
    }*/

}

void MainWindow::on_pushButtonLoad_clicked()
{
     bool flag_cancel=false;
    //QString file_name = ui->lineEditSaveDir->text();
      QString file_name = QString::null;
    // file_name += ((operating_system == win32) ? "/script.bat" : "script.sh");
    QFileInfo file_info(file_name);
    if (!file_info.exists()) {
       if(file_name_transimit==NULL)
       {
        //  qDebug("file_name");
           file_name = QFileDialog::getOpenFileName(this,
                                                    "Save as..." ,
                                                    (operating_system == os_win32 ? "script.bat" : "script.sh"),
                                                    (operating_system == os_win32 ? "script (*.bat);;Any (*.*)" : "script (*.sh);Any (*.*)"));
       if(file_name==NULL)
       {
           flag_cancel=true;
       }
       }
       else
       {
           file_name=file_name_transimit;
           file_name_transimit=QString::null;
       }
       ui->lineEditSaveDir->setText(file_name);
    }
    std::ifstream is;
    is.open(file_name.toStdString().c_str());
    std::string line;

    if( flag_cancel!=true)
    {
    ui->checkCUCorMat->setChecked(false);
    ui->checkLp->setChecked(false);
    ui->checkLp_NodalMetrics->setChecked(false);
    ui->checkCP_NodalMetrics->setChecked(false);
    ui->checkCP->setChecked(false);
  //  ui->checkSmallWordProperty->setChecked(false);
    ui->checkL_Modularity->setChecked(false);
    ui->checkPC_CPU->setChecked(false);
    ui->checkDegree->setChecked(false);
    ui->checkCUBC->setChecked(false);
    ui->checkCUEC->setChecked(false);
    ui->checkConvertNII->setChecked(false);
    ui->checkSmallWordProperty->setChecked(false);
    }
    bool flag_unweighted=false;
    bool flag_weighted=false;
    while (std::getline(is, line)) {  //这个的意思是读完整个文件，我去，太牛了！
        //ui->lineEditSaveDir->setText(line.c_str());
       // qDebug("%s",line);
        std::vector<std::string> tokens; //在命名空间中再找一个特定的命名空间
        std::string token;
        std::istringstream line_is(line);//在循环中定义的变量在每次循环里都要经历创建和撤消的过程
        //qDebug("%s",line_is);
        while (line_is >> token)
        {
            tokens.push_back(token); //把line_is的数据输到向量tokens，若输入一直有，为真，就一直输入
            qDebug("%s",token);
        }

       // qDebug("read!");
        //std::string out;
        //for (size_t i = 0; i < tokens.size(); ++i) {
        //    out += tokens[i] + " ";
        //}
        //ui->lineEditSaveDir->setText(out.c_str());

        if (tokens.empty())
        {
            continue;       //一行一行的输入！
        // qDebug("empty!");
        }
        if(tokens[0] == "echo")
        {
            //原则 1.2.若没有调用该if语句，即判断为warning
       // qDebug("read echo !");
            if(tokens[1] =="unweightednetworks")
            {

                qDebug("read echo unweighted !");
                qDebug("ffffllllaaaagggg");
               ui->switchButton->setCurrentIndex(0);  //why???换到第三行就不行？
                flag_unweighted=true;
                flag_weighted=false;//保险一点

             //   qDebug("ffffllllaaaagggg");

            }
            if((tokens[1] =="weightednetworks"))
            {
                ui->switchButton->setCurrentIndex(1);
                flag_weighted=true;
                flag_unweighted=false;

     /*           //切换界面，传递参数,发射信号
                this->hide();// 1.切换界面
                w1.setParent(this);
                w1.show();
                //2.传递参数 用全局变量extern哈哈哈！
                file_name_transimit=file_name;
                //发射信号
                emit  mySignalMtoS();
                */
            }
            if(tokens[1] =="unweightednetworks")
            {
                flag_unweighted=true;
                flag_weighted=false;//保险一点


            }
            if(tokens[1] =="unweightednetworks"||tokens[1] =="weightednetworks")
            {
                emit mySignalMclr();
               qDebug("EXECUTION");
            }

        }
        if(flag_unweighted==true&&flag_weighted==true)
        {
            QMessageBox::information(this, "Warning", "Script error(s).Error(s) may occur.Only one sort of networks type (weighted or unweighted) can be selected by using the command 'echo' in the script", QMessageBox::Ok, QMessageBox::Ok);
            //清理参数，必须这样，要不乱了·· 这个不要求echo语句的位置！好优点！``所以必须在这统一清理参数！  
            emit mySignalMclr();
        }
        if(flag_unweighted==true&&flag_weighted==false)
        {
            qDebug("red unwei specific");
       /*     if(tokens[1]!="unweightednetworks")
            {
            string temp="/unweighted";
            string tk1=tokens[1];
            qDebug("%s",tk1);
             int templen=temp.length();
             int len=tk1.length();
             int length=len-templen;
            qDebug("%d",len);
            qDebug("%d",templen);
            qDebug("%d",length);
           string unweighted=tokens[1].substr(0,tokens[1].length()-11);
             qDebug("%s",unweighted);
            }*/
            qDebug("red unwei specific");
        if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUCorMat.exe" : "./exefiles/CUCormat")) {
            if (tokens.size() >= 7) {
                ui->checkCUCorMat->setChecked(true);
                 emit mySignalMgrayBoxCUCorMat(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->mask_threshold->setText(tokens[2].c_str());
                if(tokens[3] == "n")
                {
                ui->CheckCUCorMat_to_average->setChecked(false);
                ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(false);
                ui->ordinary->setChecked(false);
                ui->fisher->setChecked(false);
                }else if(tokens[3] == "bf")
                {
                ui->CheckCUCorMat_to_average->setChecked(true);
                ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(true);
                ui->ordinary->setChecked(false);
                ui->fisher->setChecked(true);
                }else if(tokens[3] == "bn")
                {
                ui->CheckCUCorMat_to_average->setChecked(true);
                ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(true);
                ui->ordinary->setChecked(true);
                ui->fisher->setChecked(false);
                }
                ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(tokens[4] == "y");
                ui->comboBoxCUCorMat_threshold_type->setCurrentIndex(tokens[5] == "r");
                std::string token6;
                for (size_t i = 6; i < tokens.size(); ++i)
                    token6 += tokens[i] + " ";
                ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setText(token6.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe" : "./exefiles/CUBFW_Lp")) {
           string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 4) {
                if(tokens[3] =="g") {
                ui->checkLp->setChecked(true);
             //   emit mySignalMgrayBoxLp(true);
             /*   if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe" : "./exefiles/CUBFW_Lp"))
                  ui->comboBoxLp_type_for_Lp->setCurrentIndex(0);
                else if(tokens[0] ==(operating_system == os_win32 ? ".\\exefiles\\BFS_MulCPU.exe" : "./exefiles/BFS_MulCPU"))
                    ui->comboBoxLp_type_for_Lp->setCurrentIndex(1);
               */
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
                ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
                                    }

            else if(tokens[3] =="n") {
                ui->checkLp_NodalMetrics->setChecked(true);
                /*
                if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe" : "./exefiles/CUBFW_Lp"))
                  ui->comboBoxLp_type_for_Lp_NodalMetrics->setCurrentIndex(0);
                else if(tokens[0] ==(operating_system == os_win32 ? ".\\exefiles\\BFS_MulCPU.exe" : "./exefiles/BFS_MulCPU"))
                    ui->comboBoxLp_type_for_Lp_NodalMetrics->setCurrentIndex(1);
                */
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
           //     ui->lineEditLp_num_of_random_networks_NodalMetrics->setText(tokens[2].c_str());
                                    }
            else if(tokens[3] =="b") {

                ui->checkLp->setChecked(true);
            //    emit mySignalMgrayBoxLp(true);
                ui->checkLp_NodalMetrics->setChecked(true);
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
                ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
             //     ui->lineEditLp_num_of_random_networks_NodalMetrics->setText(tokens[2].c_str());
                                    }
                                    }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Cp.exe" : "./Cp")) {
            string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 4) {
                if(tokens[3]=="g"){
                ui->checkCP->setChecked(true);
                emit mySignalMgrayBoxCP(true);
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
                ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());}
            else if(tokens[3]=="n"){
                    ui->checkCP_NodalMetrics->setChecked(true);

                    ui->lineEdit_Working_Directory->setText(unweighted.c_str());
               //     ui->lineEditCp_num_of_random_networks_NodalMetrics->setText(tokens[2].c_str());
                }
                else if(tokens[3]=="b")
                {
                    ui->checkCP->setChecked(true);
                    emit mySignalMgrayBoxCP(true);
                    ui->checkCP_NodalMetrics->setChecked(true);
                    ui->lineEdit_Working_Directory->setText(unweighted.c_str());
                    ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
                }
            }
        } else if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\SmallWorldProperty.exe" : "./exefiles/SmallWorldProperty")) {
             string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 2) {
                     ui->checkSmallWordProperty->setChecked(true);
                     emit mySignalMgrayBoxSmallWorld(true);
                     ui->lineEdit_Working_Directory->setText(unweighted.c_str());
                 }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Louvain_Modularity.exe" : "./exefiles/Louvain_Modularity")
               ||tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Newman_Modularity_GPU.exe" : "./exefiles/Newman_Modularity_GPU")
               || tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Newman_Modularity_CPU.exe" : "./exefiles/Newman_Modularity_CPU")
                ) {
             string unweighted=tokens[1].substr(0,tokens[1].length()-11);
         if (tokens.size() == 3) {
             ui->checkL_Modularity->setChecked(true);
             emit mySignalMgrayBoxL_Modularity(true);
             if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Louvain_Modularity.exe" : "./exefiles/Louvain_Modularity"))
             ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(0);
             else if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Newman_Modularity_GPU.exe" : "./exefiles/Newman_Modularity_GPU") )
             ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(1);
             else ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(2);
             ui->lineEdit_Working_Directory->setText(unweighted.c_str());
             ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
         }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\PC_CPU.exe" : "./exefiles/PC_CPU")) {
             string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 2) {
                ui->checkPC_CPU->setChecked(true);
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\Degree.exe" : "./exefiles/Degree")) {
            string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 2) {
                ui->checkDegree->setChecked(true);
                emit mySignalMgrayBoxDegree(true);
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUBC.exe" : "./exefiles/CUBC")) {
             string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 2) {
                ui->checkCUBC->setChecked(true);
                emit mySignalMgrayBoxCUBC(true);
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUEC.exe" : "./exefiles/CUEC")) {
             string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 2) {
                ui->checkCUEC->setChecked(true);
                emit mySignalMgrayBoxCUEC(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\ConvertNII.exe" : "./exefiles/ConvertNII")) {
             string unweighted=tokens[1].substr(0,tokens[1].length()-11);
            if (tokens.size() == 4) {
                ui->checkConvertNII->setChecked(true);
                emit mySignalMgrayBoxConvertNII(true);
                ui->lineEdit_Working_Directory->setText(unweighted.c_str());
                ui->lineEdit_Mask_File->setText(tokens[2].c_str());
                ui->mask_threshold->setText(tokens[3].c_str());
            }
        }
        }
        if(flag_weighted==true&&flag_unweighted==false)
        {    
        if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe" : "./exefiles_weighted/CUCormat")) {
            if (tokens.size() >= 7) {
                ui->checkCUCorMat->setChecked(true);
                 emit mySignalMgrayBoxCUCorMat(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->mask_threshold->setText(tokens[2].c_str());
                if(tokens[3] == "n")
                {
                ui->CheckCUCorMat_to_average->setChecked(false);
                ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(false);
                ui->ordinary->setChecked(false);
                ui->fisher->setChecked(false);
                }else if(tokens[3] == "bf")
                {
                ui->CheckCUCorMat_to_average->setChecked(true);
                ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(true);
                ui->ordinary->setChecked(false);
                ui->fisher->setChecked(true);
                }else if(tokens[3] == "bn")
                {
                ui->CheckCUCorMat_to_average->setChecked(true);
                ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(true);
                ui->ordinary->setChecked(true);
                ui->fisher->setChecked(false);
                }
                ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(tokens[4] == "y");
                ui->comboBoxCUCorMat_threshold_type->setCurrentIndex(tokens[5] == "r");
                std::string token6;
                for (size_t i = 6; i < tokens.size(); ++i)
                    token6 += tokens[i] + " ";
                ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setText(token6.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe" : "./exefiles_weighted/CUBFW_Lp")) {
            string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 4) {
                if(tokens[3] =="g") {
                ui->checkLp->setChecked(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
                ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
                                    }

            else if(tokens[3] =="n") {
                ui->checkLp_NodalMetrics->setChecked(true);
                /*
                if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles\\CUBFW_Lp.exe" : "./exefiles/CUBFW_Lp"))
                  ui->comboBoxLp_type_for_Lp_NodalMetrics->setCurrentIndex(0);
                else if(tokens[0] ==(operating_system == os_win32 ? ".\\exefiles\\BFS_MulCPU.exe" : "./exefiles/BFS_MulCPU"))
                    ui->comboBoxLp_type_for_Lp_NodalMetrics->setCurrentIndex(1);
                */
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
           //     ui->lineEditLp_num_of_random_networks_NodalMetrics->setText(tokens[2].c_str());
                                    }
            else if(tokens[3] =="b") {

                ui->checkLp->setChecked(true);
            //    emit mySignalMgrayBoxLp(true);
                ui->checkLp_NodalMetrics->setChecked(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
                ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
             //     ui->lineEditLp_num_of_random_networks_NodalMetrics->setText(tokens[2].c_str());
            }
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe" : "./exefiles_weighted/Cp")) {
            string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 5) {
                if(tokens[4]=="g"){
                ui->checkCP->setChecked(true);
                emit mySignalMgrayBoxCP(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
                ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
                ui->comboBoxCp_Cp_type->setCurrentIndex(tokens[3] == "2");
                }
            else if(tokens[4]=="n"){
                ui->checkCP_NodalMetrics->setChecked(true);
                ui->comboBoxCp_Cp_type->setCurrentIndex(tokens[3] == "2");
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
                }
                else if(tokens[4]=="b")
                {
                    ui->checkCP->setChecked(true);
                    emit mySignalMgrayBoxCP(true);
                    ui->checkCP_NodalMetrics->setChecked(true);
                    ui->lineEdit_Working_Directory->setText(weighted.c_str());
                    ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
                    ui->comboBoxCp_Cp_type->setCurrentIndex(tokens[3] == "2");
                }
            }
        } else if(tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\SmallWorldProperty.exe" : "./exefiles_weighted/SmallWorldProperty")) {
             string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 2) {
                     ui->checkSmallWordProperty->setChecked(true);
                     ui->lineEdit_Working_Directory->setText(weighted.c_str());
                 }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\Louvain_Modularity.exe" : "./exefiles_weighted/Louvain_Modularity")
                ) {
             string weighted=tokens[1].substr(0,tokens[1].length()-9);
         if (tokens.size() == 3) {
             ui->checkL_Modularity->setChecked(true);
             emit mySignalMgrayBoxL_Modularity(true);
             ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(0);
             ui->lineEdit_Working_Directory->setText(weighted.c_str());
             ui->lineEditLp_num_of_random_networks->setText(tokens[2].c_str());
         }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\PC_CPU.exe" : "./exefiles_weighted/PC_CPU")) {
             string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 2) {
                ui->checkPC_CPU->setChecked(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\Degree.exe" : "./exefiles_weighted/Degree")) {
            string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 2) {
              //   string weighted=tokens[1].substr(0,tokens[1].length()-9);
                ui->checkDegree->setChecked(true);
                emit mySignalMgrayBoxDegree(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUEC.exe" : "./exefiles_weighted/CUEC")) {
           string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 2) {
                ui->checkCUEC->setChecked(true);
                emit mySignalMgrayBoxCUEC(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\ConvertNII.exe" : "./exefiles_weighted/ConvertNII")) {
            string weighted=tokens[1].substr(0,tokens[1].length()-9);
            if (tokens.size() == 4) {
                ui->checkConvertNII->setChecked(true);
                emit mySignalMgrayBoxConvertNII(true);
                ui->lineEdit_Working_Directory->setText(weighted.c_str());
                ui->lineEdit_Mask_File->setText(tokens[2].c_str());
                ui->mask_threshold->setText(tokens[3].c_str());
            }
        }
        }

    }

        if(flag_unweighted==false&&flag_weighted==false&&flag_cancel==false){  //这儿有debug！！！
         QMessageBox::information(this, "Warning", "Script error(s).Error(s) may occur.Please advise on the networks type (weighted or unweighted) by using the command 'echo' in the script", QMessageBox::Ok, QMessageBox::Ok);
         //清理参数，必须这样，要不乱了·· 这个不要求echo语句的位置！好优点！``所以必须在这统一清理参数！
         emit mySignalMclr();
    };

    is.close();
}


void MainWindow::on_toolButtonSaveDir_clicked()
{
    ui->lineEditSaveDir->setText(QFileDialog::getOpenFileName(this,
                                                              "Open script..." ,
                                                              (operating_system == os_win32 ? "*.bat" : "*.sh"),
                                                              (operating_system == os_win32 ? "script (*.bat);;Any (*.*)" : "script (*.sh);Any (*.*)")));
}
/*
void MainWindow::on_switchButton_clicked()
{
    this->hide();
   w1.setParent(this);
    w1.show();
}*/
void MainWindow::on_checkCUCorMat_clicked(bool checked)
{
    if(checked==0){
      //   ui->lineEditCUCorMat_Dir_for_BOLD->setEnabled(false);
      //   ui->lineEditCUCorMat_threshold_for_mask->setEnabled(false);
        // if(ui->CheckCUCorMat_to_average->isChecked())
        // ui->CheckCUCorMat_to_average->setChecked(false);
         ui->labelCUCorMat_threshold_for_correlation_coefficient->setEnabled(false);
         ui->labelCUCorMat_threshold_type->setEnabled(false);
         ui->CheckCUCorMat_to_average->setEnabled(false);
         ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(false);
         ui->ordinary->setEnabled(false);
         ui->fisher->setEnabled(false);
         ui->radioButtonCUCorMat_to_save_cormatrix->setEnabled(false);
         ui->comboBoxCUCorMat_threshold_type->setEnabled(false);
         ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setEnabled(false);
     //    ui->toolButtonCUCorMat_Dir_for_BOLD->setEnabled(false);

       //  ui->labelCUCorMat_Dir_for_BOLD->setVisible(false);
        // ui->lineEditCUCorMat_Dir_for_BOLD->setVisible(false);
       //  ui->toolButtonCUCorMat_Dir_for_BOLD->setVisible(false);
/*
         ui->labelCUCorMat_threshold_for_mask->setVisible(false);
         ui->lineEditCUCorMat_threshold_for_mask->setVisible(false);

         ui->labelCUCorMat_to_average->setVisible(false);
         ui->groupBoxCUCorMat_to_average->setVisible(false);

         //ui->labelCUCorMat_to_save_cormatrix->setVisible(false);
         ui->radioButtonCUCorMat_to_save_cormatrix->setVisible(false);

         ui->labelCUCorMat_threshold_type->setVisible(false);
         ui->comboBoxCUCorMat_threshold_type->setVisible(false);

         ui->labelCUCorMat_threshold_for_correlation_coefficient->setVisible(false);
         ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setVisible(false);
         ui->groupBox_CUCorMat->setMaximumHeight(50);
         ui->groupBox_CUCorMat->setStyleSheet("border:none");//注意！
*/
    }
     else if(checked==1)
     {
      //  ui->lineEditCUCorMat_Dir_for_BOLD->setEnabled(true);
    //    ui->lineEditCUCorMat_threshold_for_mask->setEnabled(true);
       // ui->CheckCUCorMat_to_average->setChecked(true);
        ui->labelCUCorMat_threshold_for_correlation_coefficient->setEnabled(true);
        ui->labelCUCorMat_threshold_type->setEnabled(true);
        ui->CheckCUCorMat_to_average->setEnabled(true);
        ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(false);
        ui->ordinary->setEnabled(true);
        ui->fisher->setEnabled(true);
        ui->radioButtonCUCorMat_to_save_cormatrix->setEnabled(true);
        ui->comboBoxCUCorMat_threshold_type->setEnabled(true);
        ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setEnabled(true);
     //   ui->toolButtonCUCorMat_Dir_for_BOLD->setEnabled(true);


    //    ui->labelCUCorMat_Dir_for_BOLD->setVisible(true);
    //    ui->lineEditCUCorMat_Dir_for_BOLD->setVisible(true);
    //    ui->toolButtonCUCorMat_Dir_for_BOLD->setVisible(true);
/*
        ui->labelCUCorMat_threshold_for_mask->setVisible(true);
        ui->lineEditCUCorMat_threshold_for_mask->setVisible(true);

        ui->labelCUCorMat_to_average->setVisible(true);
        ui->groupBoxCUCorMat_to_average->setVisible(true);

        //ui->labelCUCorMat_to_save_cormatrix->setVisible(true);
        ui->radioButtonCUCorMat_to_save_cormatrix->setVisible(true);

        ui->labelCUCorMat_threshold_type->setVisible(true);
        ui->comboBoxCUCorMat_threshold_type->setVisible(true);

        ui->labelCUCorMat_threshold_for_correlation_coefficient->setVisible(true);
        ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setVisible(true);
        ui->groupBox_CUCorMat->setMaximumHeight(16777215);
*/
   //     ui->groupBox_CUCorMat->setStyleSheet("QGroupBox#groupBox_CUCorMat{border: 2px solid rgb(200, 197, 191);}");//注意！

    }
}
void MainWindow::on_checkCP_clicked(bool checked)
{
    if(checked==0&&ui->switchButton->currentIndex()==1){
        ui->labelCp_Cp_type->setEnabled(false);
        ui->comboBoxCp_Cp_type->setEnabled(false);
       //  ui->lineEditCp_input_dir->setEnabled(false);
       //  ui->lineEditCp_num_of_random_networks->setEnabled(false);
        // ui->toolButtonCp_input_dir->setEnabled(false);

       //  ui->labelCp_input_dir->setVisible(false);
         //ui->lineEditCp_input_dir->setVisible(false);
      //   ui->toolButtonCp_input_dir->setVisible(false);
      //   ui->labelCp_num_of_random_networks->setVisible(false);
      //   ui->lineEditCp_num_of_random_networks->setVisible(false);
    }
     else if(checked==1&&ui->switchButton->currentIndex()==1)
     {
        ui->labelCp_Cp_type->setEnabled(true);
        ui->comboBoxCp_Cp_type->setEnabled(true);
        //ui->lineEditCp_input_dir->setEnabled(true);
       // ui->lineEditCp_num_of_random_networks->setEnabled(true);
       // ui->toolButtonCp_input_dir->setEnabled(true);

      //  ui->labelCp_input_dir->setVisible(true);
        //ui->lineEditCp_input_dir->setVisible(true);
       // ui->toolButtonCp_input_dir->setVisible(true);
      //  ui->labelCp_num_of_random_networks->setVisible(true);
      //  ui->lineEditCp_num_of_random_networks->setVisible(true);
    }
}
void MainWindow::on_checkL_Modularity_clicked(bool checked)
{
    if(checked==0&&ui->switchButton->currentIndex()==0){
         ui->comboBoxLouvain_Modularity_modularity_type->setEnabled(false);
        ui->labelL_Modularity_modularity_type->setEnabled(false);
         // ui->lineEditL_Modularity_dir_for_csr->setEnabled(false);
       //  ui->lineEditL_Modularity_num_of_random_networks->setEnabled(false);
        // ui->toolButtonL_Modularity_dir_for_csr->setEnabled(false);

      //   ui->labelL_Modularity_modularity_type->setVisible(false);
      //   ui->comboBoxLouvain_Modularity_modularity_type->setVisible(false);
        // ui->labelL_Modularity_dir_for_csr->setVisible(false);
        // ui->lineEditL_Modularity_dir_for_csr->setVisible(false);
        // ui->toolButtonL_Modularity_dir_for_csr->setVisible(false);
      //   ui->labelL_Modularity_num_of_random_networks->setVisible(false);
     //    ui->lineEditL_Modularity_num_of_random_networks->setVisible(false);
    }
     else if(checked==1&&ui->switchButton->currentIndex()==0)
     {
        ui->comboBoxLouvain_Modularity_modularity_type->setEnabled(true);
        ui->labelL_Modularity_modularity_type->setEnabled(true);
     //   ui->lineEditL_Modularity_dir_for_csr->setEnabled(true);
       // ui->lineEditL_Modularity_num_of_random_networks->setEnabled(true);
      //  ui->toolButtonL_Modularity_dir_for_csr->setEnabled(true);

      //  ui->labelL_Modularity_modularity_type->setVisible(true);
      //  ui->comboBoxLouvain_Modularity_modularity_type->setVisible(true);
      //  ui->labelL_Modularity_dir_for_csr->setVisible(true);
      //  ui->lineEditL_Modularity_dir_for_csr->setVisible(true);
      //  ui->toolButtonL_Modularity_dir_for_csr->setVisible(true);
      //  ui->labelL_Modularity_num_of_random_networks->setVisible(true);
      //  ui->lineEditL_Modularity_num_of_random_networks->setVisible(true);
    }

}

void MainWindow::on_toolButton_Working_Directory_clicked()
{
   QString directory=QFileDialog::getExistingDirectory(this, "Directory");
    ui->lineEdit_Working_Directory->setText(directory);
   // ui->lineEdit_Working_Directory->setText(directory);
    //ui->lineEditBFS_MulCPU_input_dir->setText(directory);
 //   ui->lineEditCp_input_dir->setText(directory);
  //  ui->lineEditDegree_input_dir->setText(directory);
 //   ui->lineEditCUBC_input_dir->setText(directory);
 //   ui->lineEditCUEC_input_dir->setText(directory);
 //   ui->lineEditPC_CPU_input_dir->setText(directory);
 //   ui->lineEditConvertNII_input_file->setText(directory);
}

void MainWindow::on_toolButton_Mask_File_clicked()
{
  QString file = QFileDialog::getOpenFileName(this,"NII File","","NII (*.nii)");
    ui->lineEdit_Mask_File->setText(file);
 //   ui->lineEditConvertNII_mask_file->setText(file);
}
/*
void MainWindow::on_radioButton_CUCorMat_to_average_clicked(bool checked)
{
    if(checked==0)
        individual_metrics=0;
    else if(checked==1)
        individual_metrics=1;
}
*/
void MainWindow::on_ordinary_clicked(bool checked)
{
    if(checked==0)
        average_ordinary=0;
    else if(checked==1)
    {
        average_ordinary=1;
        average_fisher=0;
    }
}

void MainWindow::on_fisher_clicked(bool checked)
{
    if(checked==0)
        average_fisher=0;
    else if(checked==1)
    {
        average_fisher=1;
        average_ordinary=0;
}
}



void MainWindow::on_radioButtonCUCorMat_to_save_cormatrix_clicked(bool checked)
{
    if(checked==0) {
        to_save_cormatrix="n";
       // ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(false);
    }
    else if(checked==1)
        to_save_cormatrix="y";
}

void MainWindow::on_CheckCUCorMat_to_average_clicked(bool checked)
{
    if(checked==0)
    {
        both_metrics=0;
    ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(false);
    ui->ordinary->setChecked(true);
    ui->fisher->setChecked(false);
    ui->ordinary->setEnabled(false);
    ui->fisher->setEnabled(false);
    }
    else if(checked==1)
    {   both_metrics=1;
        ui->groupBoxCUCorMat_to_average_groupBox->setEnabled(true);
        ui->ordinary->setChecked(true);
        ui->fisher->setChecked(false);
        ui->ordinary->setEnabled(true);
        ui->fisher->setEnabled(true);
    }
}

void MainWindow::on_switchButton_currentIndexChanged(int index)
{
    if(index==0)
    {
        ui->comboBoxCp_Cp_type->setCurrentIndex(0);
        ui->labelCp_Cp_type->setEnabled(false);
        ui->comboBoxCp_Cp_type->setEnabled(false);
        ui->labelCp_Cp_type->setToolTip("The selection of clustering coefficient type can only be activated in weighted networks");
        ui->comboBoxCp_Cp_type->setToolTip("The selection of clustering coefficient type can only be activated in weighted networks");

        ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(0);
        if(ui->checkL_Modularity->isChecked())
        {
            ui->labelL_Modularity_modularity_type->setEnabled(true);
            ui->comboBoxLouvain_Modularity_modularity_type->setEnabled(true);
        }
        ui->labelL_Modularity_modularity_type->setToolTip("We provide three specific methods for module parcellation in unweighted networks.\n"
                                                          "In algorithm ways,Louvain or Newman method could be selected. Specifically, \n"
                                                          "user can choose the implement on CPU or GPU for Newman method.");
        ui->comboBoxLouvain_Modularity_modularity_type->setToolTip("We provide three specific methods for module parcellation in unweighted networks.\n"
                                                          "In algorithm ways,Louvain or Newman method could be selected. Specifically, \n"
                                                          "user can choose the implement on CPU or GPU for Newman method.");
        ui->checkCUBC->setEnabled(true);
        ui->checkCUBC->setToolTip("The betweenness centrality is enabled in unweighted networks");
    }
    if(index==1)
    {
        ui->comboBoxCp_Cp_type->setCurrentIndex(0);
        if(ui->checkCP->isChecked())
    {
            ui->labelCp_Cp_type->setEnabled(true);
            ui->comboBoxCp_Cp_type->setEnabled(true);
    }
        ui->labelCp_Cp_type->setToolTip("In weighted networks, User can decide clustering coefficient type with different formula definition");
        ui->comboBoxCp_Cp_type->setToolTip("In weighted networks, User can decide clustering coefficient type with different formula definition");

        ui->comboBoxLouvain_Modularity_modularity_type->setCurrentIndex(0);
        ui->labelL_Modularity_modularity_type->setEnabled(false);
        ui->comboBoxLouvain_Modularity_modularity_type->setEnabled(false);

        ui->labelL_Modularity_modularity_type->setToolTip("In weighted networks");
        ui->comboBoxLouvain_Modularity_modularity_type->setToolTip("We provide three specific methods for module parcellation in unweighted networks.\n"
                                                          "In algorithm ways,Louvain or Newman method could be selected. Specifically, \n"
                                                          "user can choose the implement on CPU or GPU for Newman method.");


    ui->checkCUBC->setChecked(false);
    ui->checkCUBC->setEnabled(false);
    ui->checkCUBC->setToolTip("This algorithm can currently be used in unweighted networks");
    }
}

/*
void MainWindow::on_checkSmallWordProperty_clicked(bool checked)
{
    if(checked==0)
    {
        ui->checkSmallWordProperty_gamma->setChecked(false);
        ui->checkSmallWordProperty_lambda->setChecked(false);
        ui->checkSmallWordProperty_gamma->setEnabled(false);
        ui->checkSmallWordProperty_lambda->setEnabled(false);
    }
    if(checked==1)
    {
        ui->checkSmallWordProperty_gamma->setEnabled(true);
        ui->checkSmallWordProperty_lambda->setEnabled(true);
       
    }
}
*/
