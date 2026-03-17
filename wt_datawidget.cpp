/*
 * 文件名: wt_datawidget.cpp
 * 文件作用: 数据编辑器主窗口实现文件
 * 功能描述:
 * 1. 管理 QTabWidget，支持多标签页显示多份数据文件。
 * 2. 实现了多文件同时打开、解析与展示的功能。
 * 3. 实现了数据的同步保存与从工程文件中恢复的功能。
 * 4. 实现了 getAllDataModels，遍历所有页签收集数据模型，向下游提供数据源。
 * 5. 统一数据界面弹窗的按钮样式为“灰底黑字”，提高可视性。
 * 6. 实现了 onFilterSample 函数，集成优化的 DataSamplingDialog 弹窗进行多阶段滤波与抽样。
 * 7. 数据保存增加原生 Excel (.xlsx, .xls) 支持，利用 QXlsx 将全表完美导出。
 * 8. [重构] 彻底整合列属性定义与单位标准化，引入统一的 onUnitManager 槽函数。
 * 9. [新增] 支持用户对表格数据的直接内存操作（替换或追加），亦保留生成新文件的途径。
 */

#include "wt_datawidget.h"
#include "ui_wt_datawidget.h"
#include "modelparameter.h"
#include "dataimportdialog.h"
#include "datasamplingdialog.h"
#include "dataunitdialog.h"   // 引入全新的统一属性与单位管理对话框
#include "dataunitmanager.h"  // 引入全局物理量与单位换算系统

// 引入 QXlsx 支持数据直接写入 Excel
#include "xlsxdocument.h"
#include "xlsxformat.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QTextStream>
#include <QDateTime>

/**
 * @brief 统一应用数据对话框的样式
 * @param dialog 目标对话框指针
 * @details 规定了白底黑字、按钮灰底黑字的基础样式，符合 main.cpp 中规定的整体基调
 */
static void applyDataDialogStyle(QWidget* dialog) {
    if (!dialog) return;
    QString qss = "QWidget { color: black; background-color: white; font-family: 'Microsoft YaHei'; }"
                  "QPushButton { "
                  "   background-color: #f0f0f0; "
                  "   color: black; "
                  "   border: 1px solid #bfbfbf; "
                  "   border-radius: 3px; "
                  "   padding: 5px 15px; "
                  "   min-width: 70px; "
                  "}"
                  "QPushButton:hover { background-color: #e0e0e0; }"
                  "QPushButton:pressed { background-color: #d0d0d0; }";
    dialog->setStyleSheet(qss);
}

/**
 * @brief 构造函数
 * @param parent 父窗口指针
 */
WT_DataWidget::WT_DataWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WT_DataWidget)
{
    ui->setupUi(this);
    initUI();
    setupConnections();
}

/**
 * @brief 析构函数
 */
WT_DataWidget::~WT_DataWidget()
{
    delete ui;
}

/**
 * @brief 初始化界面
 * @details 主要用于软件启动时更新各个按钮的可用状态
 */
void WT_DataWidget::initUI()
{
    updateButtonsState();
}

/**
 * @brief 绑定界面控件的信号与槽
 * @details 集中管理所有按钮点击事件、标签页切换事件等
 */
void WT_DataWidget::setupConnections()
{
    // 文件操作区
    connect(ui->btnOpenFile, &QPushButton::clicked, this, &WT_DataWidget::onOpenFile);
    connect(ui->btnSave, &QPushButton::clicked, this, &WT_DataWidget::onSave);
    connect(ui->btnExport, &QPushButton::clicked, this, &WT_DataWidget::onExportExcel);

    // 数据处理区
    connect(ui->btnUnitManager, &QPushButton::clicked, this, &WT_DataWidget::onUnitManager); // 统一的单位管理按钮
    connect(ui->btnTimeConvert, &QPushButton::clicked, this, &WT_DataWidget::onTimeConvert);
    connect(ui->btnPressureDropCalc, &QPushButton::clicked, this, &WT_DataWidget::onPressureDropCalc);
    connect(ui->btnCalcPwf, &QPushButton::clicked, this, &WT_DataWidget::onCalcPwf);

    // 滤波与采样
    connect(ui->btnFilterSample, &QPushButton::clicked, this, &WT_DataWidget::onFilterSample);

    // 异常检查
    connect(ui->btnErrorCheck, &QPushButton::clicked, this, &WT_DataWidget::onHighlightErrors);

    // 标签页交互
    connect(ui->tabWidget, &QTabWidget::currentChanged, this, &WT_DataWidget::onTabChanged);
    connect(ui->tabWidget, &QTabWidget::tabCloseRequested, this, &WT_DataWidget::onTabCloseRequested);
}

/**
 * @brief 获取当前正在显示的单页数据表对象
 * @return 当前活动标签页的 DataSingleSheet 指针
 */
DataSingleSheet* WT_DataWidget::currentSheet() const {
    return qobject_cast<DataSingleSheet*>(ui->tabWidget->currentWidget());
}

/**
 * @brief 获取当前活动表格的数据模型
 * @return 当前的 QStandardItemModel 指针
 */
QStandardItemModel* WT_DataWidget::getDataModel() const {
    if (auto sheet = currentSheet()) return sheet->getDataModel();
    return nullptr;
}

/**
 * @brief 获取所有已打开页签的数据模型
 * @return 以文件名为键，数据模型为值的 QMap 集合，供下游模块调用
 */
QMap<QString, QStandardItemModel*> WT_DataWidget::getAllDataModels() const
{
    QMap<QString, QStandardItemModel*> map;
    for (int i = 0; i < ui->tabWidget->count(); ++i) {
        DataSingleSheet* sheet = qobject_cast<DataSingleSheet*>(ui->tabWidget->widget(i));
        if (sheet) {
            QString key = sheet->getFilePath();
            if (key.isEmpty()) key = ui->tabWidget->tabText(i);
            map.insert(key, sheet->getDataModel());
        }
    }
    return map;
}

/**
 * @brief 获取当前活动表格的文件名/路径
 */
QString WT_DataWidget::getCurrentFileName() const {
    if (auto sheet = currentSheet()) return sheet->getFilePath();
    return QString();
}

/**
 * @brief 判断当前是否已加载数据
 */
bool WT_DataWidget::hasData() const {
    return ui->tabWidget->count() > 0;
}

/**
 * @brief 更新界面按钮的启用/禁用状态
 * @details 当没有打开任何数据表时，禁用数据处理相关按钮以防止程序崩溃
 */
void WT_DataWidget::updateButtonsState()
{
    bool hasSheet = (ui->tabWidget->count() > 0);
    ui->btnSave->setEnabled(hasSheet);
    ui->btnExport->setEnabled(hasSheet);
    ui->btnUnitManager->setEnabled(hasSheet); // 同步开启/关闭单位管理按钮
    ui->btnTimeConvert->setEnabled(hasSheet);
    ui->btnPressureDropCalc->setEnabled(hasSheet);
    ui->btnCalcPwf->setEnabled(hasSheet);
    ui->btnFilterSample->setEnabled(hasSheet);
    ui->btnErrorCheck->setEnabled(hasSheet);

    if (auto sheet = currentSheet()) ui->filePathLabel->setText(sheet->getFilePath());
    else ui->filePathLabel->setText("未加载文件");
}

/**
 * @brief 打开文件槽函数
 * @details 调出文件选择框，支持 CSV、TXT、Excel 等格式的多选导入
 */
void WT_DataWidget::onOpenFile()
{
    QString filter = "所有支持文件 (*.csv *.txt *.xlsx *.xls);;Excel (*.xlsx *.xls);;CSV 文件 (*.csv);;文本文件 (*.txt);;所有文件 (*.*)";
    QStringList paths = QFileDialog::getOpenFileNames(this, "打开数据文件", "", filter);
    if (paths.isEmpty()) return;

    for (const QString& path : paths) {
        if (path.endsWith(".json", Qt::CaseInsensitive)) {
            loadData(path, "json");
            return;
        }

        DataImportDialog dlg(path, this);
        applyDataDialogStyle(&dlg);

        if (dlg.exec() == QDialog::Accepted) {
            DataImportSettings settings = dlg.getSettings();
            createNewTab(path, settings);
        }
    }
}

/**
 * @brief 创建新的数据标签页
 * @param filePath 数据文件路径
 * @param settings 导入设置（包含分隔符、编码等）
 */
void WT_DataWidget::createNewTab(const QString& filePath, const DataImportSettings& settings) {
    DataSingleSheet* sheet = new DataSingleSheet(this);
    if (sheet->loadData(filePath, settings)) {
        QFileInfo fi(filePath);
        ui->tabWidget->addTab(sheet, fi.fileName());
        ui->tabWidget->setCurrentWidget(sheet);
        connect(sheet, &DataSingleSheet::dataChanged, this, &WT_DataWidget::onSheetDataChanged);
        updateButtonsState();
        emit fileChanged(filePath, "text");
        emit dataChanged();
    } else {
        delete sheet;
        ui->statusLabel->setText("加载文件失败: " + filePath);
    }
}

/**
 * @brief 根据路径和类型加载数据
 * @param filePath 文件路径
 * @param fileType 文件类型
 */
void WT_DataWidget::loadData(const QString& filePath, const QString& fileType)
{
    if (fileType == "json") return; // 工程文件另行处理
    DataImportDialog dlg(filePath, this);
    applyDataDialogStyle(&dlg);
    if (dlg.exec() == QDialog::Accepted) {
        createNewTab(filePath, dlg.getSettings());
    }
}

/**
 * @brief 保存所有数据到工程文件
 * @details 遍历所有标签页转化为 JSON 格式，并调用全局模型参数实例进行工程保存
 */
void WT_DataWidget::onSave() {
    QJsonArray allData;
    for (int i = 0; i < ui->tabWidget->count(); ++i) {
        DataSingleSheet* sheet = qobject_cast<DataSingleSheet*>(ui->tabWidget->widget(i));
        if (sheet) allData.append(sheet->saveToJson());
    }
    ModelParameter::instance()->saveTableData(allData);
    ModelParameter::instance()->saveProject();

    QMessageBox msgBox(this);
    msgBox.setWindowTitle("保存");
    msgBox.setText("所有标签页数据已同步保存到项目文件。");
    msgBox.setIcon(QMessageBox::Information);
    msgBox.addButton(QMessageBox::Ok);
    applyDataDialogStyle(&msgBox);
    msgBox.exec();
}

/**
 * @brief 从工程数据恢复表格视图
 * @details 当打开工程文件时，调用此函数将 JSON 还原为多个标签页数据
 */
void WT_DataWidget::loadFromProjectData() {
    clearAllData();
    QJsonArray dataArray = ModelParameter::instance()->getTableData();
    if (dataArray.isEmpty()) {
        ui->statusLabel->setText("无数据");
        return;
    }

    bool isNewFormat = false;
    if (!dataArray.isEmpty()) {
        QJsonValue first = dataArray.first();
        if (first.isObject() && first.toObject().contains("filePath") && first.toObject().contains("data")) {
            isNewFormat = true;
        }
    }

    if (isNewFormat) {
        for (auto val : dataArray) {
            QJsonObject sheetObj = val.toObject();
            DataSingleSheet* sheet = new DataSingleSheet(this);
            sheet->loadFromJson(sheetObj);
            QString path = sheet->getFilePath();
            QFileInfo fi(path);
            ui->tabWidget->addTab(sheet, fi.fileName().isEmpty() ? "恢复数据" : fi.fileName());
            connect(sheet, &DataSingleSheet::dataChanged, this, &WT_DataWidget::onSheetDataChanged);
        }
    } else {
        DataSingleSheet* sheet = new DataSingleSheet(this);
        QJsonObject sheetObj;
        sheetObj["filePath"] = "Restored Data";
        if (!dataArray.isEmpty() && dataArray.first().toObject().contains("headers")) {
            sheetObj["headers"] = dataArray.first().toObject()["headers"];
        }
        QJsonArray rows;
        for (int i = 1; i < dataArray.size(); ++i) {
            QJsonObject rowObj = dataArray[i].toObject();
            if (rowObj.contains("row_data")) rows.append(rowObj["row_data"]);
        }
        sheetObj["data"] = rows;
        sheet->loadFromJson(sheetObj);
        ui->tabWidget->addTab(sheet, "恢复数据");
        connect(sheet, &DataSingleSheet::dataChanged, this, &WT_DataWidget::onSheetDataChanged);
    }

    updateButtonsState();
    ui->statusLabel->setText("数据已恢复");
}

/**
 * @brief 清空所有已加载的数据与标签页
 */
void WT_DataWidget::clearAllData() {
    ui->tabWidget->clear();
    ui->filePathLabel->setText("未加载文件");
    ui->statusLabel->setText("无数据");
    updateButtonsState();
    emit dataChanged();
}

/**
 * @brief 核心槽函数：统一单位管理与列属性定义
 * @details 通过 DataUnitDialog 的配置结果，分流处理为“直接内存修改当前表”或“另存为新文件”两套机制。
 */
void WT_DataWidget::onUnitManager()
{
    QStandardItemModel* model = getDataModel();
    if (!model || model->rowCount() == 0) {
        QMessageBox::warning(this, "提示", "当前无数据可操作！");
        return;
    }

    DataUnitDialog dialog(model, this);
    if (dialog.exec() != QDialog::Accepted) return;

    auto tasks = dialog.getConversionTasks();
    if (tasks.isEmpty()) return;

    bool appendMode = dialog.isAppendMode();
    bool saveToFile = dialog.isSaveToFile();
    int rowCount = model->rowCount();
    int colCount = model->columnCount();

    if (!saveToFile) {
        // ========== 模式 A：直接修改当前内存数据表格（不生成文件） ==========
        if (appendMode) {
            // 追加模式：在模型最后插入新列
            int insertPos = colCount;
            for (const auto& task : tasks) {
                model->insertColumn(insertPos);
                model->setHeaderData(insertPos, Qt::Horizontal, task.newHeaderName);
                for (int r = 0; r < rowCount; ++r) {
                    QString origText = model->item(r, task.colIndex)->text();
                    if (task.needsMathConversion) {
                        bool ok; double val = origText.toDouble(&ok);
                        if (ok) {
                            double newVal = DataUnitManager::instance()->convert(val, task.quantityType, task.fromUnit, task.toUnit);
                            model->setItem(r, insertPos, new QStandardItem(QString::number(newVal, 'g', 6)));
                        } else {
                            model->setItem(r, insertPos, new QStandardItem(""));
                        }
                    } else {
                        model->setItem(r, insertPos, new QStandardItem(origText)); // 仅定义属性时直接复制数值
                    }
                }
                insertPos++;
            }
        } else {
            // 替换模式：覆盖原列的数据
            for (const auto& task : tasks) {
                model->setHeaderData(task.colIndex, Qt::Horizontal, task.newHeaderName);
                if (task.needsMathConversion) {
                    for (int r = 0; r < rowCount; ++r) {
                        QString origText = model->item(r, task.colIndex)->text();
                        bool ok; double val = origText.toDouble(&ok);
                        if (ok) {
                            double newVal = DataUnitManager::instance()->convert(val, task.quantityType, task.fromUnit, task.toUnit);
                            model->item(r, task.colIndex)->setText(QString::number(newVal, 'g', 6));
                        }
                    }
                }
            }
        }
        emit dataChanged(); // 发射更新信号供绘图及其他核心模块刷新
        ui->statusLabel->setText("列属性与单位已在当前表格中更新完毕。");

    } else {
        // ========== 模式 B：生成新文件并加载为新页签 ==========
        QVector<QStringList> fullTable;
        QStringList headers;

        // 收集表头
        for (int c = 0; c < colCount; ++c) {
            headers.append(model->horizontalHeaderItem(c) ? model->horizontalHeaderItem(c)->text() : "");
        }

        if (appendMode) {
            for (const auto& task : tasks) headers.append(task.newHeaderName);
        } else {
            for (const auto& task : tasks) headers[task.colIndex] = task.newHeaderName;
        }

        // 遍历并换算数值内容
        for (int r = 0; r < rowCount; ++r) {
            QStringList rowData;
            for (int c = 0; c < colCount; ++c) rowData.append(model->item(r, c)->text());

            if (appendMode) {
                for (const auto& task : tasks) {
                    if (task.needsMathConversion) {
                        bool ok; double val = rowData[task.colIndex].toDouble(&ok);
                        if (ok) {
                            double nV = DataUnitManager::instance()->convert(val, task.quantityType, task.fromUnit, task.toUnit);
                            rowData.append(QString::number(nV, 'g', 6));
                        } else {
                            rowData.append("");
                        }
                    } else {
                        rowData.append(rowData[task.colIndex]); // 照常追加原始值
                    }
                }
            } else {
                for (const auto& task : tasks) {
                    if (task.needsMathConversion) {
                        bool ok; double val = rowData[task.colIndex].toDouble(&ok);
                        if (ok) {
                            double nV = DataUnitManager::instance()->convert(val, task.quantityType, task.fromUnit, task.toUnit);
                            rowData[task.colIndex] = QString::number(nV, 'g', 6);
                        }
                    }
                }
            }
            fullTable.append(rowData);
        }

        // 调用复用的稳健落地函数
        QString currentPath = getCurrentFileName();
        saveAndLoadNewData(currentPath, headers, fullTable);
        ui->statusLabel->setText("列属性与单位已应用，并生成新数据文件。");
    }
}

/**
 * @brief 滤波与取样核心更新槽函数
 * @details 唤出滤波窗口，处理完成后通过 saveAndLoadNewData 重新挂载数据
 */
void WT_DataWidget::onFilterSample()
{
    QStandardItemModel* model = getDataModel();
    if (!model || model->rowCount() == 0) {
        QMessageBox::warning(this, "提示", "当前无数据可处理！");
        return;
    }

    DataSamplingDialog dialog(model, this);
    if (dialog.exec() == QDialog::Accepted) {

        QVector<QStringList> processedTable = dialog.getProcessedTable();
        QStringList headers = dialog.getHeaders();

        if (processedTable.isEmpty()) {
            QMessageBox::warning(this, "错误", "处理后数据为空，请检查参数设置！");
            return;
        }

        QString currentPath = getCurrentFileName();
        saveAndLoadNewData(currentPath, headers, processedTable);
    }
}

/**
 * @brief 辅助函数：将内存中的数据表写入本地文件并重新加载为新标签页
 * @param oldFilePath 原始文件路径，用于提取目录和前缀
 * @param headers 表头字符串列表
 * @param fullTable 数据主体的二维字符串列表
 * @details 支持保存为 Excel (.xlsx) / CSV / TXT 格式
 */
void WT_DataWidget::saveAndLoadNewData(const QString& oldFilePath, const QStringList& headers, const QVector<QStringList>& fullTable)
{
    QFileInfo fi(oldFilePath);
    QString baseName = fi.completeBaseName();
    QString dir = fi.absolutePath();

    // 生成默认文件名
    QString defaultName = QString("%1_处理后_%2.xlsx").arg(baseName).arg(QDateTime::currentDateTime().toString("HHmmss"));
    QString defaultPath = dir + "/" + defaultName;

    // 弹窗增加了 Excel 选项
    QString savePath = QFileDialog::getSaveFileName(this, "保存处理后的全表数据", defaultPath, "Excel文件 (*.xlsx *.xls);;CSV文件 (*.csv);;文本文件 (*.txt)");
    if (savePath.isEmpty()) return;

    bool isExcel = savePath.endsWith(".xlsx", Qt::CaseInsensitive) || savePath.endsWith(".xls", Qt::CaseInsensitive);

    if (isExcel) {
        // ========== 写入 Excel 文件 ==========
        QXlsx::Document xlsx;

        // 1. 写入头文件信息
        QXlsx::Format infoFormat;
        infoFormat.setFontColor(Qt::darkGray);
        xlsx.write(1, 1, "// PWT System: Data Processed/Standardized Data", infoFormat);
        xlsx.write(2, 1, "// Original File: " + fi.fileName(), infoFormat);

        // 2. 写入表头
        QXlsx::Format headerFormat;
        headerFormat.setFontBold(true);
        headerFormat.setHorizontalAlignment(QXlsx::Format::AlignHCenter);
        headerFormat.setPatternBackgroundColor(QColor(240, 240, 240));

        for (int c = 0; c < headers.size(); ++c) {
            xlsx.write(3, c + 1, headers[c], headerFormat);
        }

        // 3. 写入数据内容
        for (int r = 0; r < fullTable.size(); ++r) {
            for (int c = 0; c < fullTable[r].size(); ++c) {
                bool ok;
                double val = fullTable[r][c].toDouble(&ok);
                if (ok) xlsx.write(r + 4, c + 1, val);
                else xlsx.write(r + 4, c + 1, fullTable[r][c]); // 非数字写文本
            }
        }

        if (!xlsx.saveAs(savePath)) {
            QMessageBox::critical(this, "错误", "无法保存 Excel 文件，请确保文件未被其它程序占用。");
            return;
        }
    } else {
        // ========== 写入 CSV/TXT 文件 ==========
        QFile file(savePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QMessageBox::critical(this, "错误", "无法创建或写入文件！");
            return;
        }
        QTextStream out(&file);
        out << "// PWT System: Data Processed/Standardized Data\n";
        out << "// Original File: " << fi.fileName() << "\n";
        out << headers.join(",") << "\n";
        for (const QStringList& row : fullTable) {
            out << row.join(",") << "\n";
        }
        file.close();
    }

    // 重新导入生成的新文件
    DataImportSettings settings;
    settings.filePath = savePath;
    settings.encoding = "UTF-8";
    settings.separator = ",";
    settings.useHeader = true;
    settings.headerRow = 3;
    settings.startRow = 4;
    settings.isExcel = isExcel; // 标记是否为 Excel 供解析模块使用

    createNewTab(savePath, settings);

    ui->statusLabel->setText("处理完成，保留有效点数: " + QString::number(fullTable.size()));
}

// ============ 单页操作委托调用区域 ============
// 以下操作全部向下穿透委托给当前正在显示的 DataSingleSheet 页面执行

/**
 * @brief 导出为 Excel
 */
void WT_DataWidget::onExportExcel() { if (auto s = currentSheet()) s->onExportExcel(); }

/**
 * @brief 执行时间格式转换计算
 */
void WT_DataWidget::onTimeConvert() { if (auto s = currentSheet()) s->onTimeConvert(); }

/**
 * @brief 计算压降
 */
void WT_DataWidget::onPressureDropCalc() { if (auto s = currentSheet()) s->onPressureDropCalc(); }

/**
 * @brief 计算井底流压 Pwf
 */
void WT_DataWidget::onCalcPwf() { if (auto s = currentSheet()) s->onCalcPwf(); }

/**
 * @brief 高亮标出表格内的异常数据（如非数字、空行等）
 */
void WT_DataWidget::onHighlightErrors() { if (auto s = currentSheet()) s->onHighlightErrors(); }

// ============ 标签页及信号联动机制 ============

/**
 * @brief 当用户切换不同的数据标签页时触发
 * @param index 切换后的页签索引
 */
void WT_DataWidget::onTabChanged(int index) {
    Q_UNUSED(index);
    updateButtonsState(); // 动态更新顶部工具栏按钮的可点击状态
    emit dataChanged();   // 通知下游模块当前显示的数据源可能发生变更
}

/**
 * @brief 当用户点击标签页右上角关闭按钮时触发
 * @param index 欲关闭的页签索引
 */
void WT_DataWidget::onTabCloseRequested(int index) {
    QWidget* widget = ui->tabWidget->widget(index);
    if (widget) {
        ui->tabWidget->removeTab(index);
        delete widget;
    }
    updateButtonsState();
    emit dataChanged();
}

/**
 * @brief 捕捉从 DataSingleSheet 发送上来的数据修改信号
 */
void WT_DataWidget::onSheetDataChanged() {
    // 只有当信号发送者是当前显示的页面时，才向外扩散数据变更信号
    if (sender() == currentSheet()) {
        emit dataChanged();
    }
}
