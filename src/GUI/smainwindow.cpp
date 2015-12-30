#include "smainwindow.h"
#include "ui_smainwindow.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QtGlobal>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
using  std::string;

QString file_name_transimit=QString::null ;


bool s_individual_metrics=0;
bool s_average_ordinary=0;
bool s_average_fisher=0;
string s_to_save_cormatrix="n";

sMainWindow::sMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::sMainWindow)
{
    ui->setupUi(this);
      pWid = NULL;
   //   this->setWindowFlags(windowFlags()& ~Qt::WindowMaximizeButtonHint);

      connect(this,SIGNAL(mySignalSclr()),this,SLOT(clearscreen()));

      connect(this,SIGNAL(mySignalSgrayBoxCUCorMat(bool)),this,SLOT(on_checkCUCorMat_clicked(bool)));
      connect(this,SIGNAL(mySignalSgrayBoxCUBFW_Lp(bool)),this,SLOT(on_checkCUBFW_Lp_clicked(bool)));
      //connect(this,SIGNAL(mySignalSgrayBoxCUEC(bool)),this,SLOT(on_checkCUEC_clicked(bool)));
      connect(this,SIGNAL(mySignalSgrayBoxCP(bool)),this,SLOT(on_checkCP_clicked(bool)));
      //connect(this,SIGNAL(mySignalSgrayBoxDegree(bool)),this,SLOT(on_checkDegree_clicked(bool)));
      //connect(this,SIGNAL(mySignalSgrayBoxCUCP(bool)),this,SLOT(on_checkCUCP_clicked(bool)));
      connect(this,SIGNAL(mySignalSgrayBoxL_Modularity(bool)),this,SLOT(on_checkL_Modularity_clicked(bool)));
      connect(this,SIGNAL(mySignalSgrayBoxPC_CPU(bool)),this,SLOT(on_checkPC_CPU_clicked(bool)));
      connect(this,SIGNAL(mySignalSgrayBoxConvertNII(bool)),this,SLOT(on_checkConvertNII_clicked(bool)));
      connect(this,SIGNAL(mySignalSgrayBoxCP_Nodal_Metrics(bool)),this,SLOT(on_checkCP_Nodal_Metrics_clicked(bool)));

     emit mySignalSclr();
    ui->statusbar->hide();
}

sMainWindow::~sMainWindow()
{
    delete ui;
}

void sMainWindow::on_pushButtonSave_clicked()
{
    if (operating_system != os_win32 && operating_system != os_linux) {
        std::cout << "Not win32 or linux!" << std::endl;

    }

    std::stringstream script;
    script <<  "echo weightednetworks" << std::endl;
    if (ui->checkCUCorMat->isChecked()) {
        //组合之外的统统选n
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditCUCorMat_threshold_for_mask->text().isEmpty()
                //|| ui->lineEditCUCorMat_to_average->text().isEmpty()
                //|| ui->lineEditCUCorMat_s_to_save_cormatrix->text().isEmpty()
                //|| ui->lineEditCUCorMat_threshold_type->text().isEmpty()
                || ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string  type=ui->comboBoxCUCorMat_threshold_type->currentText().toStdString()=="correlation"?"r":"s";
        if(s_individual_metrics==0&&s_average_ordinary==1&&s_average_fisher==0)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_mask->text().toStdString() <<
              ' ' <<
              "yn" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              s_to_save_cormatrix <<
                  //ui->lineEditCUCorMat_s_to_save_cormatrix->text().toStdString() <<
              ' ' <<
              type<<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else if(s_individual_metrics==1&&s_average_ordinary==1&&s_average_fisher==0)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_mask->text().toStdString() <<
              ' ' <<
              "bn" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              s_to_save_cormatrix <<
                  //ui->lineEditCUCorMat_s_to_save_cormatrix->text().toStdString() <<
              ' ' <<
              type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else if(s_individual_metrics==0&&s_average_ordinary==1&&s_average_fisher==1)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_mask->text().toStdString() <<
              ' ' <<
              "yf" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              s_to_save_cormatrix <<
                  //ui->lineEditCUCorMat_s_to_save_cormatrix->text().toStdString() <<
              ' ' <<
             type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else if(s_individual_metrics==1&&s_average_ordinary==1&&s_average_fisher==1)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_mask->text().toStdString() <<
              ' ' <<
              "bf" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              s_to_save_cormatrix <<
                  //ui->lineEditCUCorMat_s_to_save_cormatrix->text().toStdString() <<
              ' ' <<
             type <<

              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
        else //if(s_individual_metrics==1&&s_average_ordinary==1&&s_average_fisher==0)
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe " : "./exefiles_weighted/CUCormat ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_mask->text().toStdString() <<
              ' ' <<
              "n" <<
                  //ui->lineEditCUCorMat_to_average->text().toStdString() <<
              ' ' <<
              s_to_save_cormatrix <<
                  //ui->lineEditCUCorMat_s_to_save_cormatrix->text().toStdString() <<
              ' ' <<
              type <<
                  //ui->lineEditCUCorMat_threshold_type->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCorMat_threshold_for_correlation_coefficient->text().toStdString() <<
              std::endl;
      // delete(type.c_str());
    }

    if (ui->checkCUBFW_Lp->isChecked()&&(!ui->checkCUBFW_Lp_Nodal_Metrics->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditCUBFW_Lp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe " : "./exefiles_weighted/CUBFW_Lp ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUBFW_Lp_num_of_random_networks->text().toStdString() <<' '<<"g"<<
              std::endl;
    }



  /*  if (ui->checkCUCP->isChecked()) {
        if (ui->lineEditCUCP_input_dir->text().isEmpty()
                || ui->lineEditCUCP_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCP.exe " : "./exefiles_weighted/CUCP ") <<
              ui->lineEditCUCP_input_dir->text().toStdString() <<
              ' ' <<
              ui->lineEditCUCP_num_of_random_networks->text().toStdString() <<
              std::endl;
    }*/

    if (ui->checkCP->isChecked()&&(!ui->checkCP_Nodal_Metrics->isChecked())) {
        if (ui->lineEditCp_num_of_random_networks->text().isEmpty()
                || ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
         string s2=ui->comboBoxCp_Cp_type->currentText().toStdString();
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCp_num_of_random_networks->text().toStdString() <<' ' <<(s2 == "Onnela" ? "2 " : " 1 ")<<' '<<"g"<<
              std::endl;
    }

    if (ui->checkDegree->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Degree.exe " : "./exefiles_weighted/Degree ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              std::endl;
    }

    if (ui->checkCUEC->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUEC.exe " : "./exefiles_weighted/CUBC ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              std::endl;
    }



    if (ui->checkL_Modularity->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditL_Modularity_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }                                           //这句话是看是不是win32 是就是一种格式 不是就是另一种格式  \\是c++中的转义字符  代表win中的反斜杠；系统路径，编程路径 \\或/ 理解：因为\会被编译错

     //   string s=ui->comboBoxLouvain_Modularity_modularity_type->currentText().toStdString();
     //   if(s=="Louvain")
           script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Louvain_Modularity.exe " : "./exefiles_weighted/Louvain_Modularity ") <<
                 ui->lineEdit_Working_Directory->text().toStdString() <<
                 ' ' <<
                 ui->lineEditL_Modularity_num_of_random_networks->text().toStdString() <<
                 std::endl;
   /*     else if (s=="Newman_GPU")
            script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Newman_Modularity_GPU.exe " : "./exefiles_weighted/Newman_Modularity_GPU ") <<
                        ui->lineEditL_Modularity_dir_for_csr->text().toStdString() <<
                        ' ' <<
                        ui->lineEditL_Modularity_num_of_random_networks->text().toStdString() <<
                        std::endl;
         else
            script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Newman_Modularity_CPU.exe " : "./exefiles_weighted/Newman_Modularity_CPU ") <<
                         ui->lineEditL_Modularity_dir_for_csr->text().toStdString() <<
                         ' ' <<
                         ui->lineEditL_Modularity_num_of_random_networks->text().toStdString() <<
                         std::endl;
   */

    }
    if (ui->checkPC_CPU->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
             //   || ui->lineEditPC_CPU_type_for_participant_coefficient->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        string s1=ui->comboBoxPC_CPU_type_for_participant_coefficient->currentText().toStdString();
        if(s1=="normalized")
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\PC_CPU.exe " : "./exefiles_weighted/PC_CPU ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
           //   ui->lineEditPC_CPU_type_for_participant_coefficient->text().toStdString() <<
              "n" <<
                  std::endl;
        else  script << (operating_system == os_win32 ? ".\\exefiles_weighted\\PC_CPU.exe " : "./exefiles_weighted/PC_CPU ") <<
                        ui->lineEdit_Working_Directory->text().toStdString() <<
                        ' ' <<
                            std::endl;
    }
    if (ui->checkConvertNII->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEdit_Mask_File->text().isEmpty()
                || ui->lineEditConvertNII_mask_threshold->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\ConvertNII.exe " : "./exefiles_weighted/ConvertNII ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEdit_Mask_File->text().toStdString() <<
              ' ' <<
              ui->lineEditConvertNII_mask_threshold->text().toStdString() <<
              std::endl;
    }
    if (ui->checkCUBFW_Lp_Nodal_Metrics->isChecked()&&(!ui->checkCUBFW_Lp->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                ) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe " : "./exefiles_weighted/CUBFW_Lp ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              "0" <<' '<<"n"<<
              std::endl;
    }

    if (ui->checkCUBFW_Lp->isChecked()&&ui->checkCUBFW_Lp_Nodal_Metrics->isChecked()) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()
                || ui->lineEditCUBFW_Lp_num_of_random_networks->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe " : "./exefiles_weighted/CUBFW_Lp ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCUBFW_Lp_num_of_random_networks->text().toStdString() <<' '<<"b"<<
              std::endl;
    }
    if (ui->checkCP_Nodal_Metrics->isChecked()&&(!ui->checkCP->isChecked())) {
        if (ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
         string s2=ui->comboBoxCp_Cp_type_Nodal_Metrics->currentText().toStdString();
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              "0"<<' ' <<(s2 == "Onnela" ? "2 " : " 1 ")<<' '<<"n"<<
              std::endl;
    }
    if (ui->checkCP_Nodal_Metrics->isChecked()&&ui->checkCP->isChecked()) {
        if (ui->lineEditCp_num_of_random_networks->text().isEmpty()
                ||ui->lineEdit_Working_Directory->text().isEmpty()) {
            QMessageBox::information(this, "Error", "Empty parameter(s).", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        if(ui->comboBoxCp_Cp_type->currentText().toStdString()!=ui->comboBoxCp_Cp_type_Nodal_Metrics->currentText().toStdString()) {
            QMessageBox::information(this, "Error", "If Clustering Coefficient is selected in global and nodal metrics at the same time, you must ensure that the 'type' parameters are the same", QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
         string s2=ui->comboBoxCp_Cp_type->currentText().toStdString();
        script << (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe " : "./exefiles_weighted/Cp ") <<
              ui->lineEdit_Working_Directory->text().toStdString() <<
              ' ' <<
              ui->lineEditCp_num_of_random_networks->text().toStdString()
               <<' ' <<(s2 == "Onnela" ? "2 " : " 1 ")<<' '<<"b"<<
              std::endl;
    }
    QString file_name = ui->lineEditSaveDir->text();
    std::cout << script.str() << std::endl;
    QFileInfo file_info(file_name);
    if (!file_info.exists()) {
        file_name = QFileDialog::getSaveFileName(this,
                                                         "Save as..." ,
                                                         (operating_system == os_win32 ? "script.bat" : "script.sh"),
                                                         (operating_system == os_win32 ? "script (*.bat);;Any (*.*)" : "script (*.sh);Any (*.*)"));

    }
    if (!file_name.isNull()) {
        ui->lineEditSaveDir->setText(file_name);
        std::ofstream os;
        os.open(file_name.toStdString().c_str());
        os << script.str();
        os.close();

        if (operating_system == os_win32) {
            std::string cmd = "start " + file_name.toStdString() + " &";
            system(cmd.c_str());
        } else {
            std::string cmd = "sh " + file_name.toStdString() + " &";
            system(cmd.c_str());
        }
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


void sMainWindow::on_pushButtonLoad_clicked()
{
    bool flag_cancel=false;
    //QString file_name = ui->lineEditSaveDir->text();
    QString file_name = QString::null ;
    // file_name += ((operating_system == win32) ? "/script.bat" : "script.sh");
    QFileInfo file_info(file_name);
    if (!file_info.exists()) {
        if(file_name_transimit==NULL){
        file_name = QFileDialog::getOpenFileName(this,
                                                 "Save as..." ,
                                                 (operating_system == os_win32 ? "script.bat" : "script.sh"),
                                                 (operating_system == os_win32 ? "script (*.bat);;Any (*.*)" : "script (*.sh);Any (*.*)"));
        if(file_name==NULL){
            flag_cancel=true;
        }
        }else {
            file_name=file_name_transimit;
            file_name_transimit=QString::null ;
            qDebug("finish assignment in childWindow");
        }
         ui->lineEditSaveDir->setText(file_name);
    }
    std::ifstream is;
    is.open(file_name.toStdString().c_str());
    std::string line;

    ui->checkCUCorMat->setChecked(false);
    ui->checkCUBFW_Lp->setChecked(false);
    ui->checkPC_CPU->setChecked(false);
   // ui->checkCUCP->setChecked(false);
    ui->checkCP->setChecked(false);
    ui->checkDegree->setChecked(false);
    ui->checkCUEC->setChecked(false);
    ui->checkConvertNII->setChecked(false);
    ui->checkL_Modularity->setChecked(false);
    int flag=0;
    while (std::getline(is, line)) {
        //ui->lineEditSaveDir->setText(line.c_str());
        std::vector<std::string> tokens; //在命名空间中再找一个特定的命名空间
        std::string token;
        std::istringstream line_is(line);
        while (line_is >> token) {
            tokens.push_back(token);
        }

        //std::string out;
        //for (size_t i = 0; i < tokens.size(); ++i) {
        //    out += tokens[i] + " ";
        //}
        //ui->lineEditSaveDir->setText(out.c_str());

        if (tokens.empty())
            continue;
        if(tokens[0] =="echo"){  //原则 1.2.若没有调用该if语句，即判断为warning

            if(tokens[1] =="weightednetworks"){
                flag++;
                //不切换界面，处理参数
                continue;
            }else if((tokens[1] =="unweightednetworks")){
                flag++;
                //切换界面，传递参数
                // 1.切换界面
                this->close();
                    if(pWid)
                    {
                        pWid->showNormal();
                        qDebug("event close");
                    }
                //2.传递参数 用全局变量extern哈哈哈！
                file_name_transimit=file_name;
                //3.发射信号
                emit  mySignalStoM();
            }
        }

        if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCorMat.exe" : "./exefiles_weighted/CUCormat")) {
            if (tokens.size() >= 7) {
                ui->checkCUCorMat->setChecked(true);
                 emit mySignalSgrayBoxCUCorMat(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->lineEditCUCorMat_threshold_for_mask->setText(tokens[2].c_str());
                if(tokens[3] == "n")
                {
                ui->radioButton_CUCorMat_to_average->setChecked(true);
                ui->ordinary->setChecked(false);
                ui->fisher->setChecked(false);
                }else if(tokens[3] == "yf")
                {
                ui->radioButton_CUCorMat_to_average->setChecked(false);
                ui->ordinary->setChecked(true);
                ui->fisher->setChecked(true);
                }else if(tokens[3] == "yn")
                {
                ui->radioButton_CUCorMat_to_average->setChecked(false);
                ui->ordinary->setChecked(true);
                ui->fisher->setChecked(false);
                }else if(tokens[3] == "bn")
                {
                ui->radioButton_CUCorMat_to_average->setChecked(true);
                ui->ordinary->setChecked(true);
                ui->fisher->setChecked(false);
                }else if(tokens[3] == "bf")
                {
                ui->radioButton_CUCorMat_to_average->setChecked(true);
                ui->ordinary->setChecked(true);
                ui->fisher->setChecked(true);
                }
                ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(tokens[4] == "y");
                ui->comboBoxCUCorMat_threshold_type->setCurrentIndex(tokens[5] == "r");
                //ui->lineEditCUCorMat_to_average->setText(tokens[3].c_str());
                //ui->lineEditCUCorMat_to_save_cormatrix->setText(tokens[4].c_str());
                //ui->lineEditCUCorMat_threshold_type->setText(tokens[5].c_str());
                std::string token6;
                for (size_t i = 6; i < tokens.size(); ++i)
                    token6 += tokens[i] + " ";
                ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setText(token6.c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUBFW_Lp.exe" : "./exefiles_weighted/CUBFW_Lp")) {
            if (tokens.size() == 4) {
              if(tokens[3] == "g") {
                ui->checkCUBFW_Lp->setChecked(true);
                emit mySignalSgrayBoxCUBFW_Lp(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->lineEditCUBFW_Lp_num_of_random_networks->setText(tokens[2].c_str());
              }
              else if(tokens[3] == "n") {
                ui->checkCUBFW_Lp_Nodal_Metrics->setChecked(true);
                //emit mySignalSgrayBoxCUBFW_Lp(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                //ui->lineEditCUBFW_Lp_num_of_random_networks->setText(tokens[2].c_str());
              }
              else if(tokens[3] == "b") {
                  ui->checkCUBFW_Lp->setChecked(true);
                  emit mySignalSgrayBoxCUBFW_Lp(true);
                ui->checkCUBFW_Lp_Nodal_Metrics->setChecked(true);
                //emit mySignalSgrayBoxCUBFW_Lp(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->lineEditCUBFW_Lp_num_of_random_networks->setText(tokens[2].c_str());
              }
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\PC_CPU.exe" : "./exefiles_weighted/PC_CPU")) {
            if (tokens.size() == 3) {
                ui->checkPC_CPU->setChecked(true);
                emit mySignalSgrayBoxPC_CPU(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->comboBoxPC_CPU_type_for_participant_coefficient->setCurrentIndex(tokens[2] == "n");
            }
            else if (tokens.size() == 2) {
                ui->checkPC_CPU->setChecked(true);
                emit mySignalSgrayBoxPC_CPU(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->comboBoxPC_CPU_type_for_participant_coefficient->setCurrentIndex(0);
            }
        } /*else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUCP.exe" : "./exefiles_weighted/CUCP")) {
            if (tokens.size() == 3) {
                ui->checkCUCP->setChecked(true);
                emit mySignalSgrayBoxCUCP(true);
                ui->lineEditCUCP_input_dir->setText(tokens[1].c_str());
                ui->lineEditCUCP_num_of_random_networks->setText(tokens[2].c_str());
            }
        }*/ else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\Cp.exe" : "./Cp")) {
            if (tokens.size() == 5) {
              if(tokens[4] == "g") {
                ui->checkCP->setChecked(true);
                emit mySignalSgrayBoxCP(true);
              //  qDebug("enter cp");
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->lineEditCp_num_of_random_networks->setText(tokens[2].c_str());
                ui->comboBoxCp_Cp_type->setCurrentIndex(tokens[3] == "1");
              }
              else if(tokens[4] == "n") {
              ui->checkCP_Nodal_Metrics->setChecked(true);
              emit mySignalSgrayBoxCP_Nodal_Metrics(true);
            //  qDebug("enter cp");
              ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
           //   ui->lineEditCp_num_of_random_networks->setText(tokens[2].c_str());
              ui->comboBoxCp_Cp_type_Nodal_Metrics->setCurrentIndex(tokens[3] == "1");
              }
              else if(tokens[4] == "b") {
                  ui->checkCP->setChecked(true);
                  emit mySignalSgrayBoxCP(true);
                  ui->checkCP_Nodal_Metrics->setChecked(true);
                  emit mySignalSgrayBoxCP_Nodal_Metrics(true);
                //  qDebug("enter cp");
                  ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                  ui->lineEditCp_num_of_random_networks->setText(tokens[2].c_str());
                  ui->comboBoxCp_Cp_type->setCurrentIndex(tokens[3] == "1");
                  ui->comboBoxCp_Cp_type_Nodal_Metrics->setCurrentIndex(tokens[3] == "1");
              }
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\Degree.exe" : "./exefiles_weighted/Degree")) {
            if (tokens.size() == 2) {
                ui->checkDegree->setChecked(true);
                emit mySignalSgrayBoxDegree(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
            }

        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\CUEC.exe" : "./exefiles_weighted/CUEC")) {
            if (tokens.size() == 2) {
                ui->checkCUEC->setChecked(true);
                emit mySignalSgrayBoxCUEC(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\ConvertNII.exe" : "./exefiles_weighted/ConvertNII")) {
            if (tokens.size() == 4) {
                ui->checkConvertNII->setChecked(true);
                emit mySignalSgrayBoxConvertNII(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->lineEdit_Mask_File->setText(tokens[2].c_str());
                ui->lineEditConvertNII_mask_threshold->setText(tokens[3].c_str());
            }
        } else if (tokens[0] == (operating_system == os_win32 ? ".\\exefiles_weighted\\Louvain_Modularity.exe" : "./exefiles_weighted/Louvain_Modularity")) {
            if (tokens.size() == 3) {
                ui->checkL_Modularity->setChecked(true);
                emit mySignalSgrayBoxL_Modularity(true);
                ui->lineEdit_Working_Directory->setText(tokens[1].c_str());
                ui->lineEditL_Modularity_num_of_random_networks->setText(tokens[2].c_str());
            }
        }

        // ui->lineEditSaveDir->setText(QString(token.c_str()));
    }
    if(flag==0&&flag_cancel==false){
         QMessageBox::information(this, "Warning", "Empty parameter(s).Error(s) may occur.Please advise on the network type (weighted or unweighted) by using the command:echo", QMessageBox::Ok, QMessageBox::Ok);
         //清理参数，必须这样，要不乱了·· 这个不要求echo语句的位置！好优点！``所以必须在这统一清理参数！
         emit mySignalSclr();
    };
    is.close();

}


/*
void sMainWindow::on_toolButtonCUCorMat_Dir_for_BOLD_clicked()
{
 ui->lineEditCUCorMat_Dir_for_BOLD->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}
*/
/*
void sMainWindow::on_toolButtonCUBFW_Lp_input_dir_clicked()
{
 ui->lineEditCUBFW_Lp_input_dir->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}

void sMainWindow::on_toolButtonPC_CPU_input_dir_clicked()
{
     ui->lineEditPC_CPU_input_dir->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}
*//*
void sMainWindow::on_toolButtonCUCP_input_dir_clicked()
{
    ui->lineEditCUCP_input_dir->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}*/
/*
void sMainWindow::on_toolButtonCp_input_dir_clicked()
{
    ui->lineEditCp_input_dir->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}

void sMainWindow::on_toolButtonDegree_input_dir_clicked()
{
    ui->lineEditDegree_input_dir->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}

void sMainWindow::on_toolButtonCUEC_input_dir_clicked()
{
     ui->lineEditCUEC_input_dir->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}

void sMainWindow::on_toolButtonConvertNII_input_file_clicked()
{
    ui->lineEditConvertNII_input_file->setText(QFileDialog::getExistingDirectory(this,
                                                                            "Directory"
                                                                            ));
}
*/
/*
void sMainWindow::on_toolButtonConvertNII_mask_file_clicked()
{
    ui->lineEditConvertNII_mask_file->setText(QFileDialog::getExistingDirectory(this,
                                                                              "Directory"

                                                                              ));
    ui->lineEditConvertNII_mask_file->setText(QFileDialog::getOpenFileName(this,
                                                                          "NII File",
                                                                          "",
                                                                          "NII (*.nii)"));
}
*/
/*
void sMainWindow::on_toolButtonL_Modularity_dir_for_csr_clicked()
{
    ui->lineEditL_Modularity_dir_for_csr->setText(QFileDialog::getExistingDirectory(this, "Directory"));
}
*/
void sMainWindow::on_toolButtonSaveDir_clicked()
{
    ui->lineEditSaveDir->setText(QFileDialog::getOpenFileName(this,
                                                              "Open script..." ,
                                                              (operating_system == os_win32 ? "*.bat" : "*.sh"),
                                                              (operating_system == os_win32 ? "script (*.bat);;Any (*.*)" : "script (*.sh);Any (*.*)")));
}
void sMainWindow::setParent(QWidget *parent){

pWid = parent;
}
void sMainWindow::on_switchButton_clicked()
{
this->close();
    if(pWid)
    {
        pWid->showNormal();
        qDebug("event close");
    }
}



void sMainWindow::on_checkCUCorMat_clicked(bool checked)
{
    if(checked==0){
         ui->lineEditCUCorMat_threshold_for_mask->setEnabled(false);
         ui->radioButton_CUCorMat_to_average->setEnabled(false);
         ui->ordinary->setEnabled(false);
         ui->fisher->setEnabled(false);
         ui->radioButtonCUCorMat_to_save_cormatrix->setEnabled(false);
         ui->comboBoxCUCorMat_threshold_type->setEnabled(false);
         ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setEnabled(false);
    }
     else if(checked==1)
     {
      //  ui->lineEditCUCorMat_Dir_for_BOLD->setEnabled(true);
        ui->lineEditCUCorMat_threshold_for_mask->setEnabled(true);
        ui->radioButton_CUCorMat_to_average->setEnabled(true);
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
        ui->groupBox_CUCorMat->setStyleSheet("QGroupBox#groupBox_CUCorMat{border: 2px solid rgb(200, 197, 191);}");//注意！
*/       // ui->groupBox_CUCorMat->setStyleSheet("font:14");//注意
    }
}

void sMainWindow::on_checkCUBFW_Lp_clicked(bool checked)
{
    if(checked==0){
     //    ui->lineEditCUBFW_Lp_input_dir->setEnabled(false);
         ui->lineEditCUBFW_Lp_num_of_random_networks->setEnabled(false);
     //    ui->toolButtonCUBFW_Lp_input_dir->setEnabled(false);
    }
     else if(checked==1)
     {
      //  ui->lineEditCUBFW_Lp_input_dir->setEnabled(true);
        ui->lineEditCUBFW_Lp_num_of_random_networks->setEnabled(true);
      //  ui->toolButtonCUBFW_Lp_input_dir->setEnabled(true);
    }
}

/*void sMainWindow::on_checkCUCP_clicked(bool checked)
{
    if(checked==0){
         ui->lineEditCUCP_input_dir->setEnabled(false);
         ui->lineEditCUCP_num_of_random_networks->setEnabled(false);
         ui->toolButtonCUCP_input_dir->setEnabled(false);
    }
     else if(checked==1)
     {
        ui->lineEditCUCP_input_dir->setEnabled(true);
        ui->lineEditCUCP_num_of_random_networks->setEnabled(true);
        ui->toolButtonCUCP_input_dir->setEnabled(true);
    }
}*/

void sMainWindow::on_checkCP_clicked(bool checked)
{
    if(checked==0){
       //  ui->lineEditCp_input_dir->setEnabled(false);
         ui->lineEditCp_num_of_random_networks->setEnabled(false);
       //  ui->toolButtonCp_input_dir->setEnabled(false);
          ui->comboBoxCp_Cp_type->setEnabled(false);
    }
     else if(checked==1)
     {
        //ui->lineEditCp_input_dir->setEnabled(true);
        ui->lineEditCp_num_of_random_networks->setEnabled(true);
       // ui->toolButtonCp_input_dir->setEnabled(true);
        ui->comboBoxCp_Cp_type->setEnabled(true);
    }
}
/*
void sMainWindow::on_checkDegree_clicked(bool checked)
{
    if(checked==0){
         ui->lineEditDegree_input_dir->setEnabled(false);
         ui->toolButtonDegree_input_dir->setEnabled(false);
    }
     else if(checked==1)
     {
        ui->lineEditDegree_input_dir->setEnabled(true);
        ui->toolButtonDegree_input_dir->setEnabled(true);
     }
}

void sMainWindow::on_checkCUEC_clicked(bool checked)
{
    if(checked==0){
       //  ui->lineEditCUEC_input_dir->setEnabled(false);
         ui->toolButtonCUEC_input_dir->setEnabled(false);
    }
     else if(checked==1)
     {
        ui->lineEditCUEC_input_dir->setEnabled(true);
        ui->toolButtonCUEC_input_dir->setEnabled(true);
    }
}
*/
void sMainWindow::on_checkL_Modularity_clicked(bool checked)
{
    if(checked==0){
      //   ui->comboBoxLouvain_Modularity_modularity_type->setEnabled(false);
   //      ui->lineEditL_Modularity_dir_for_csr->setEnabled(false);
         ui->lineEditL_Modularity_num_of_random_networks->setEnabled(false);
     //    ui->toolButtonL_Modularity_dir_for_csr->setEnabled(false);
    }
     else if(checked==1)
     {
       // ui->comboBoxLouvain_Modularity_modularity_type->setEnabled(true);
       // ui->lineEditL_Modularity_dir_for_csr->setEnabled(true);
        ui->lineEditL_Modularity_num_of_random_networks->setEnabled(true);
    //    ui->toolButtonL_Modularity_dir_for_csr->setEnabled(true);
    }

}

void sMainWindow::on_checkPC_CPU_clicked(bool checked)
{
    if(checked==0){
      //   ui->lineEditPC_CPU_input_dir->setEnabled(false);
      //   ui->toolButtonPC_CPU_input_dir->setEnabled(false); //注意这个名称的命名两个界面不一样！！
         ui->comboBoxPC_CPU_type_for_participant_coefficient->setEnabled(false);
    }
     else if(checked==1)
     {
      //  ui->lineEditPC_CPU_input_dir->setEnabled(true);
      //  ui->toolButtonPC_CPU_input_dir->setEnabled(true);
        ui->comboBoxPC_CPU_type_for_participant_coefficient->setEnabled(true);
    }
}

void sMainWindow::on_checkConvertNII_clicked(bool checked)
{
    if(checked==0){
        // ui->lineEditConvertNII_input_file->setEnabled(false);
         //ui->lineEditConvertNII_mask_file->setEnabled(false);
         ui->lineEditConvertNII_mask_threshold->setEnabled(false);
         //ui->toolButtonConvertNII_input_file->setEnabled(false);
    //     ui->toolButtonConvertNII_mask_file->setEnabled(false);
    }
     else if(checked==1)
     {
       // ui->lineEditConvertNII_input_file->setEnabled(true);
      //  ui->lineEditConvertNII_mask_file->setEnabled(true);
        ui->lineEditConvertNII_mask_threshold->setEnabled(true);
      //  ui->toolButtonConvertNII_input_file->setEnabled(true);
        //ui->toolButtonConvertNII_mask_file->setEnabled(true);
    }

}
void sMainWindow::clearscreen(){
    //清屏函数两个功能：1.清空内容；2.使界面回到灰框状态
    //1.清空内容；
    ui->lineEditCUCorMat_threshold_for_mask->setText("");
    ui->radioButton_CUCorMat_to_average->setChecked(false);
    ui->ordinary->setChecked(false);
    ui->fisher->setChecked(false);

    ui->radioButtonCUCorMat_to_save_cormatrix->setChecked(false);
    ui->comboBoxCUCorMat_threshold_type->setCurrentIndex(0);
    ui->lineEditCUCorMat_threshold_for_correlation_coefficient->setText("");
   // ui->lineEditCUBFW_Lp_input_dir->setText("");;
    ui->lineEditCUBFW_Lp_num_of_random_networks->setText("");
   // ui->toolButtonCUBFW_Lp_input_dir->setEnabled(false);
    ui->lineEdit_Working_Directory->setText("");;
  //  ui->lineEditLp_num_of_random_networks->setText("");
    //ui->lineEditCp_input_dir->setText("");
    ui->lineEditCp_num_of_random_networks->setText("");
    ui->comboBoxCp_Cp_type->setCurrentIndex(0);
    ui->comboBoxCp_Cp_type_Nodal_Metrics->setCurrentIndex(0);

//    ui->lineEditCUCP_input_dir->setText("");
//    ui->lineEditCUCP_num_of_random_networks->setText("");

      //ui->lineEditDegree_input_dir->setText("");

       //ui->lineEditCUEC_input_dir->setText("");

       //ui->lineEditL_Modularity_dir_for_csr->setText("");
       ui->lineEditL_Modularity_num_of_random_networks->setText("");

      // ui->lineEditPC_CPU_input_dir->setText("");
       ui->comboBoxPC_CPU_type_for_participant_coefficient->setCurrentIndex(0);

       //ui->lineEditConvertNII_input_file->setText("");
    //   ui->lineEditConvertNII_mask_file->setText("");
       ui->lineEditConvertNII_mask_threshold->setText("");
       //2.使界面回到灰框状态
       emit mySignalSgrayBoxCUCorMat(false);
       emit mySignalSgrayBoxCUBFW_Lp(false);
       emit mySignalSgrayBoxCUEC(false);
       emit mySignalSgrayBoxCP(false);
       emit mySignalSgrayBoxDegree(false);
       emit mySignalSgrayBoxCUCP(false);
       emit mySignalSgrayBoxL_Modularity(false);
       emit mySignalSgrayBoxPC_CPU(false);
       emit mySignalSgrayBoxConvertNII(false);
       emit mySignalSgrayBoxCP_Nodal_Metrics(false);

}

void sMainWindow::on_toolButton_Working_Directory_clicked()
{
    QString directory=QFileDialog::getExistingDirectory(this, "Directory");
     ui->lineEdit_Working_Directory->setText(directory);
}

void sMainWindow::on_toolButton_Mask_File_clicked()
{
    QString file = QFileDialog::getOpenFileName(this,"NII File","","NII (*.nii)");
      ui->lineEdit_Mask_File->setText(file);
    //  ui->lineEditConvertNII_mask_file->setText(file);
}

void sMainWindow::on_checkCP_Nodal_Metrics_clicked(bool checked)
{
    if(checked==0){
       //  ui->lineEditCp_input_dir->setEnabled(false);
       //  ui->lineEditCp_num_of_random_networks->setEnabled(false);
       //  ui->toolButtonCp_input_dir->setEnabled(false);
          ui->comboBoxCp_Cp_type_Nodal_Metrics->setEnabled(false);
    }
     else if(checked==1)
     {
        //ui->lineEditCp_input_dir->setEnabled(true);
       // ui->lineEditCp_num_of_random_networks->setEnabled(true);
       // ui->toolButtonCp_input_dir->setEnabled(true);
        ui->comboBoxCp_Cp_type_Nodal_Metrics->setEnabled(true);
    }
}

void sMainWindow::on_radioButton_CUCorMat_to_average_clicked(bool checked)
{
    if(checked==0)
        s_individual_metrics=0;
    else if(checked==1)
        s_individual_metrics=1;
}

void sMainWindow::on_ordinary_clicked(bool checked)
{
    if(checked==0)
        s_average_ordinary=0;
    else if(checked==1)
        s_average_ordinary=1;
}

void sMainWindow::on_fisher_clicked(bool checked)
{
    if(checked==0)
        s_average_fisher=0;
    else if(checked==1)
        s_average_fisher=1;

}


void sMainWindow::on_radioButtonCUCorMat_to_save_cormatrix_clicked(bool checked)
{

    if(checked==0) {
        s_to_save_cormatrix="n";
       // ui->radioButtonCUCorMat_s_to_save_cormatrix->setChecked(false);
    }
    else if(checked==1)
        s_to_save_cormatrix="y";
}
