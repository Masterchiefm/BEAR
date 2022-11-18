import os
import json

from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication, QTableWidgetItem, QFileDialog, QMessageBox
import pandas as pd

from gui import Ui_MainWindow

from analyser import SangerBaseCall

import requests
def getLyric():
    try:
        url2 = 'https://v1.jinrishici.com/all'
        lyric = requests.get(url2, timeout=1).json()
        content = lyric['content']
        try:
            origin = lyric['origin']
        except:
            origin = "Unknown"
        try:
            author = lyric['author']
        except:
            author = "Unknown"


        output = content + "\n\t\t\t" + "——《" + origin + "》\t" + author
        return output
    except Exception as e:
        output = "BaseEdit-Analyser" + "\n\t\t\t" + "——Written by M.Q. at ShanghaiTech University"
        return output


class MyMainWin(QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None):
        super(MyMainWin, self).__init__(parent)

        self.setupUi(self)
        self.version = "1.2.0"
        self.setWindowTitle(self.windowTitle() + "v"+self.version)
        output = "BaseEdit-Analyser" + "\n\t\t\t" + "——Written by M.Q. at ShanghaiTech University"
        self.label_lyric.setText(output)
        self.selected_rows = []

        col_count = self.tableWidget.columnCount()
        self.col_names = []
        self.col_name_locations = {}
        for i in range(col_count):
            name = self.tableWidget.horizontalHeaderItem(i).text()
            self.col_names.append(name)
            self.col_name_locations[name] = i

        print(self.col_names)


        # self.tableWidget.resizeColumnToContents(0)
        self.tableWidget.selectColumn(0)

        self.tableWidget.itemSelectionChanged.connect(self.showSelection)
        self.pushButton_add_line.clicked.connect(self.addLine)
        self.pushButton_clear_table.clicked.connect(self.clearTable)

        self.pushButton_import_from_sheet.clicked.connect(self.importFromSheet)
        self.pushButton_export_sheet.clicked.connect(self.exportSheet)
        self.pushButton_start.clicked.connect(self.annalyse)

        self.plainTextEdit_sanger_rec.dropped.connect(self.inputSanger)

        self.pushButton_open_report.clicked.connect(self.openReport)

        self.pushButton_open_folder.clicked.connect(self.openFolder)

        self.pushButton_del_lines.clicked.connect(self.delLine)

        self.plainTextEdit_sanger_rec_auto.dropped.connect(self.autoInputSanger)

        self.tableWidget.clicked.connect(self.disableAutoFill)

        self.actionEnglish.triggered.connect(self.loadTanslation)
        self.actionAbout.triggered.connect(self.showAbout)
        self.actionTutor.triggered.connect(self.showTutor)
        self.checkUpdate()


    def checkUpdate(self):
        """使用requests模块和GitHub api获取最新版本"""

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'}
        latestInfo = requests.get("https://gitee.com/api/v5/repos/MasterChiefm/BEAR/releases/latest", timeout=2,  headers= headers)
        # print(latestInfo)
        try:
            info = latestInfo.text
            info = json.loads(info)
            # print(info)
            latestVersion = info["tag_name"]
            releaseInfo = info["body"]
            #print(latestVersion)

            if latestVersion == self.version:
                # QMessageBox.about(self, "更新", "已经是最新")
                return
            else:
                msg = latestVersion + "\n" + releaseInfo
                # QMessageBox.about(self,"更新",msg)
                answer = QMessageBox.question(self,"New BEAR version found", msg + "\n\n Update now?")
                if answer == QMessageBox.Yes:
                    os.startfile("https://gitee.com/MasterChiefm/BEAR/releases/latest")
                else:
                    return
        except Exception as e:
            # QMessageBox.about(self,"网络错误","无法连接到GitHub服务器")
            print(e)

    def showAbout(self):
        os.startfile("https://gitee.com/MasterChiefm/BEAR/blob/master/README.md")

    def showTutor(self):
        os.startfile("https://www.bilibili.com/video/BV1sG411E71G")

    def loadTanslation(self):
        path, type_ = QFileDialog.getOpenFileName(self, "Load translation file","","Qt language file(*.qm)")
        try:
            trans = QtCore.QTranslator()
            trans.load(path)
            _app = QApplication.instance()
            _app.installTranslator(trans)
            self.retranslateUi(self)
        except Exception as e:
            QMessageBox.about(self, "Translate UI Fail","Fail to change laguage!\n\n" + str(e))

    def toChinese(self):
        trans = QtCore.QTranslator()
        # trans.load("./en")
        _app = QApplication.instance()
        _app.installTranslator(trans)



    def disableAutoFill(self):
        self.checkBox_auto_fill_col.setChecked(False)

    def autoInputSanger(self):
        lyric = getLyric()
        self.label_lyric.setText(lyric)

        file_list = self.plainTextEdit_sanger_rec_auto.toPlainText().split("file:///")

        self.plainTextEdit_sanger_rec_auto.clear()
        separator = self.lineEdit_separator.text().strip()
        # msg = True
        for file in file_list:
            file_path = file.replace("\n","")
            if os.path.isfile(file_path):
                file_name = os.path.basename(file_path)
                if "ab1" in file_name.lower():
                    sample_name = file_name.split(separator)[0]

                    # if msg:
                    #     message = 'You are using auto fill function.\n\nThe abi file will be separated by "' +\
                    #               separator + '"\n\n' + "For example, your \n\n" +\
                    #               file_name + '\n\nwill be regarded as "' + \
                    #               sample_name + '" to fill in the table.\n\nAre you sure to do it?'
                    #     answer = QMessageBox.question(self,"Auto fill table", message)
                    #     # print(answer)
                    #     msg = False
                    #     if answer == QMessageBox.Yes:
                    #         pass
                    #     else:
                    #         return

                    exist = False
                    for row in range(self.tableWidget.rowCount()):
                        sample_exist = self.tableWidget.item(row,0)
                        try:
                            if sample_name == sample_exist.text():
                                self.tableWidget.setItem(row, self.col_name_locations["Sequencing File"],QTableWidgetItem(file_path))
                                exist = True
                        except:
                            continue

                    if exist == False:
                        newRow = self.tableWidget.rowCount() + 1
                        self.tableWidget.setRowCount(newRow)
                        self.tableWidget.setItem(newRow-1, self.col_name_locations["Sample name"], QTableWidgetItem(sample_name))
                        self.tableWidget.setItem(newRow-1, self.col_name_locations["Sequencing File"], QTableWidgetItem(file_path))



    def openFolder(self):
        os.startfile(self.lineEdit_path.text())


    def openReport(self):

        selection = self.selected_rows
        for i in selection:
            try:
                report_file = self.tableWidget.item(i, self.col_name_locations["Report"]).text()
                os.startfile(report_file)
            except Exception as e:
                QMessageBox.about(self,"ERROR","Report Not Found!\n" + str(e))


    def delLine(self):
        selection = self.selected_rows
        for i in selection:
            try:
                self.tableWidget.removeRow(selection[0])
            except Exception as e:
                print(e)




    def exportSheet(self,tem_save = False):
        lyric = getLyric()
        self.label_lyric.setText(lyric)

        if tem_save == True:
            file_path = self.lineEdit_path.text() + "/临时存储.xlsx"
        else:
            file_path, type = QFileDialog.getSaveFileName(self, "存", "", "excel(*.xlsx)")

        sheet = pd.DataFrame(columns = self.col_names)
        if file_path:
            pass
        else:
            return

        table = self.tableWidget

        for row in range(table.rowCount()):
            data = []
            for col in range(len(self.col_names)):
                try:
                    text = table.item(row, col).text()
                except:
                    text = ""
                data.append(text)
            sheet.loc[row] = data

        sheet.to_excel(file_path)


    def importFromSheet(self):
        file_path, type = QFileDialog.getOpenFileName(self,"导入", "","excel(*.xlsx)")
        if file_path:
            pass
        else:
            return

        sheet = pd.read_excel(file_path, index_col = 0)
        # count = self.tableWidget.rowCount()
        self.tableWidget.clearContents()
        self.tableWidget.setRowCount(len(sheet.index))

        for row in range(len(sheet.index)):
            data = sheet.iloc[row]
            for col in range(len(sheet.columns)):
                text = str(sheet.iloc[row][col])
                # if col == 0:
                #     if text == "nan":
                #         return
                if text == "nan":
                    text = ""

                self.tableWidget.setItem(row,col,QTableWidgetItem(text))









    def annalyse(self):

        self.label_lyric.setText("❤❤❤❤❤❤   Analyzing... Please waite!   ❤❤❤❤❤❤")
        lyric = getLyric()
        if self.setSavePath():
            pass
        else:
            self.label_lyric.setText(lyric)
            return
        table = self.tableWidget
        # sheet = pd.DataFrame(columns=self.col_names)

        for row in range(self.tableWidget.rowCount()):
            try:
                name = table.item(row, 0).text().strip()
                sg = table.item(row, 1).text().strip()
                try:
                    base_from = table.item(row,2).text().upper().strip()
                except:
                    base_from = ""
                try:
                    base_to = table.item(row,3).text().upper().strip()
                except:
                    base_to = ""

                sanger_file = table.item(row, 5).text()
                print(base_to)

                if name == "":
                    continue
                if sg == "":
                    continue
                if sanger_file == "":
                    continue

            except Exception as e:
                QMessageBox.about(self,"ERROR","please check input" + str(e))
                continue



            # 开始分析

            try:
                sanger = SangerBaseCall(sanger_file)
                result_sheet = sanger.getTargetData(target_seq=sg,annalyse_window=5)
                # if result_sheet == "None":
                #     QMessageBox.about(self, "ERROR", "Target Not Found!")



                try:
                    # 这个没有填写的话就不分析了
                    base_position = int(table.item(row, 4).text().strip())
                    perc = result_sheet.loc[base_to, base_position]
                    self.tableWidget.setItem(row,self.col_name_locations["Edit efficiency"],QTableWidgetItem(str("%.2f%%" % perc)))
                except:
                    pass
            except Exception as e:
                # QMessageBox.about(self,"ERROR",sanger_file + "error:\n" + str(e))
                self.tableWidget.setItem(row, self.col_name_locations["Edit efficiency"], QTableWidgetItem("Error, please check your input" + str(e)))
                continue
            if type(result_sheet) == type("a"):
                self.tableWidget.setItem(row, self.col_name_locations["Edit efficiency"], QTableWidgetItem("Error, target not found!"))
                continue

            for i in result_sheet.columns:
                try:
                    perc = result_sheet.loc[base_to, i]
                    if base_from == result_sheet.loc["sequence",i]:
                        self.tableWidget.setItem(row, self.col_name_locations[str(i)], QTableWidgetItem(str("%.2f%%" % perc)))
                except Exception as e:
                    print(e)

            report_path = self.lineEdit_path.text()+"/reports/" + name + ".html"
            try:
                base_to = table.item(row, 3).text().upper().strip()
            except:
                base_to = ""

            sanger.generateReport(sg, name, report_path,base_of_interest=base_to)
            self.tableWidget.setItem(row, self.col_name_locations["Report"], QTableWidgetItem(str(report_path)))



            #self.exportSheet(tem_save = True)
        msg = "Done! \nThe html reports can be found in \n" + self.lineEdit_path.text() + "/reports"
        self.label_lyric.setText(lyric)
        QMessageBox.about(self,"Done!", msg)







    def setSavePath(self):
        save_path = QFileDialog.getExistingDirectory(self,"选路径")
        if save_path:
            pass
        else:
            return False
        self.lineEdit_path.setText(save_path)
        try:
            os.mkdir(save_path + "/" + "reports")
        except:
            pass

        return True





    def addLine(self):
        current_count = self.tableWidget.rowCount()
        self.tableWidget.setRowCount(current_count + 1)
        new_row = current_count + 1
        self.tableWidget.setItem(new_row,0,QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 1, QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 2, QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 3, QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 4, QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 5, QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 6, QTableWidgetItem(""))
        self.tableWidget.setItem(new_row, 7, QTableWidgetItem(""))

    def showSelection(self):
        selection = self.tableWidget.selectedIndexes()
        # self.currentSelectedIndex = selection
        # self.tabWidget.setCurrentIndex(0)
        rows = []
        names = []

        for i in selection:
            row = i.row()
            try:
                name = self.tableWidget.item(row, 0).text()
            except:
                name = ""

            if row in rows:
                pass
            else:
                rows.append(row)
                names.append(name)
        if len(rows) > 10:
            self.label_selection.setText("选中了" + str(len(rows)) + "个")
        else:
            self.label_selection.setText(str(names))
        self.selected_rows = rows

        col = self.tableWidget.currentColumn()
        if self.checkBox_auto_fill_col.isChecked():
            try:
                first_item = self.tableWidget.item(rows[0], col)
                first_data = first_item.text()
            except:
                first_data = ""

            for i in selection:
                row = i.row()
                self.tableWidget.setItem(row, col, QTableWidgetItem(first_data))




    def clearTable(self):
        self.tableWidget.clearContents()
        self.tableWidget.setRowCount(0)
        lyric = getLyric()
        self.label_lyric.setText(lyric)


    # def clearRec(self):
    #     self.plainTextEdit_control_rec.clear()
    #     self.plainTextEdit_reference_rec.clear()
    #     self.plainTextEdit_sample_rec.clear()

    def inputSanger(self):
        text = self.plainTextEdit_sanger_rec.toPlainText()
        text = text.replace("file:///","")
        for i in self.selected_rows:
            text_item = QTableWidgetItem(text)
            self.tableWidget.setItem(i, 5, text_item)

        self.plainTextEdit_sanger_rec.clear()






if __name__ == "__main__":
    import sys

    # trans = QtCore.QTranslator()
    # trans.load("./en")


    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    # app.installTranslator(trans)
    win = MyMainWin()
    win.show()
    sys.exit(app.exec_())