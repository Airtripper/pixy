//
// begin license header
//
// This file is part of Pixy CMUcam5 or "Pixy" for short
//
// All Pixy source code is provided under the terms of the
// GNU General Public License v2 (http://www.gnu.org/licenses/gpl-2.0.html).
// Those wishing to use Pixy source code, software and/or
// technologies under different licensing terms should contact us at
// cmucam@cs.cmu.edu. Such licensing terms are available for
// all portions of the Pixy codebase presented here.
//
// end license header
//

#include "interpreter.h"
#include "configdialog.h"
#include "pixytypes.h"
#include "paramfile.h"
#include <QLabel>
#include <QMessageBox>
#include "debug.h"
#include <QTableWidget>
#include <QPushButton>
#include <QAbstractButton>
#include <QCheckBox>
#include <QSlider>
#include <QFileDialog>
#include <stdexcept>



ConfigDialog::ConfigDialog(QWidget *parent, Interpreter *interpreter) : QDialog(parent), m_ui(new Ui::ConfigDialog)
{

    m_ui->setupUi(this);
    m_interpreter = interpreter;
    m_interpreter->unwait(); // unhang interpreter if it's waiting

    m_pixyTabs = new QTabWidget(this);
    m_ui->pixyLayout->addWidget(m_pixyTabs);
    m_pixymonTabs = new QTabWidget(this);
    m_ui->pixymonLayout->addWidget(m_pixymonTabs);
    connect(interpreter, SIGNAL(paramLoaded()), this, SLOT(load()));
    connect(m_ui->buttonBox, SIGNAL(clicked(QAbstractButton*)), this, SLOT(apply(QAbstractButton*)));

    m_interpreter->loadParams();

    render(m_interpreter->m_pixymonParameters, m_ui->pixymonLayout, m_pixymonTabs);

#ifdef __MACOS__
    setMinimumWidth(650);
#endif
#ifdef __LINUX__
    setMinimumWidth(600);
#endif
#ifdef __WINDOWS__
    setMinimumWidth(550);
#endif
}


ConfigDialog::~ConfigDialog()
{
}

QWidget *ConfigDialog::findCategory(const QString &category, QTabWidget *tabs)
{
    int i;
    QString tabText;

    for (i=0; true; i++)
    {
        tabText = tabs->tabText(i);
        if (tabText=="")
            break;
        else if (tabText==category)
            return tabs->widget(i);
    }
    return NULL;
}

int ConfigDialog::updateDB()
{
    int res1, res2;

    res1 = updateDB(&m_interpreter->m_pixyParameters);
    res2 = updateDB(m_interpreter->m_pixymonParameters);

    if (res1<0)
        return res1;
    if (res2<0)
        return res2;
    return 0;
}

int ConfigDialog::updateDB(ParameterDB *data)
{
    Parameters &parameters = data->parameters();

    for (int i=0; i<parameters.size(); i++)
    {
        Parameter &parameter = parameters[i];
        int type = parameter.type();
        uint flags = parameter.property(PP_FLAGS).toInt();

        parameter.clearShadow();

        if (flags&PRM_FLAG_INTERNAL) // don't render!
            continue;

        QLineEdit *line = (QLineEdit *)parameter.property(PP_WIDGET).toLongLong();

        if (type==PT_FLT32)
        {
            bool ok;
            float val;
            val = line->text().toFloat(&ok);
            if (!ok)
            {
                QMessageBox::critical(NULL, "Error", parameter.id() + " needs to be a floating point number!");
                return -1;
            }
            if (val!=parameter.value().toFloat())
            {
                parameter.set(val);
                parameter.setDirty(true);
            }
            // handle slider if applicable
            if (flags&PRM_FLAG_SLIDER)
            {
                float min = parameter.property(PP_MIN).toFloat();
                float max = parameter.property(PP_MAX).toFloat();
                QSlider *slider = (QSlider *)parameter.property(PP_WIDGET2).toLongLong();
                float pos;
                // value = min + pos/100(max-min)
                //(value - min)100/(max-min) = pos
                pos = (val - min)*SLIDER_SIZE/(max - min);
                slider->setSliderPosition((int)pos);
            }
        }
        else if (type==PT_STRING)
        {
            QString val = line->text();
            if (flags&PRM_FLAG_PATH)
            {
                QDir dir(val);

                if (!dir.exists())
                {
                    QMessageBox::critical(NULL, "Error", parameter.id() + " \"" + val + "\" is not a valid folder!");
                    return -1;
                }
            }
            if (val!=parameter.value().toString())
            {
                parameter.set(val);
                parameter.setDirty(true);
            }
        }
        else // must be int type
        {
            if (flags&PRM_FLAG_CHECKBOX) // checkbox is a special case of int
            {
                QCheckBox *cbox = (QCheckBox *)line;

                if(cbox->isChecked()!=parameter.value().toBool())
                {
                    parameter.set(cbox->isChecked());
                    parameter.setDirty(true);
                }
            }
            else
            {

                int base;
                bool ok;
                if (line->text().left(2)=="0x")
                    base = 16;
                else
                    base = 10;
                if (flags&PRM_FLAG_SIGNED)
                {
                    int val = line->text().toInt(&ok, base);
                    if (!ok)
                    {
                        QMessageBox::critical(NULL, "Error", parameter.id() + " needs to be an integer!");
                        return -1;
                    }
                    if (val!=parameter.valueInt())
                    {
                        parameter.set(val);
                        parameter.setDirty(true);
                    }
                }
                else
                {
                    uint val = line->text().toUInt(&ok, base);
                    if (!ok)
                    {
                        QMessageBox::critical(NULL, "Error", parameter.id() + " needs to be an unsigned integer!");
                        return -1;
                    }
                    if (val!=parameter.value().toUInt())
                    {
                        parameter.set(val);
                        parameter.setDirty(true);
                    }
                }
                // handle slider if applicable
                if (flags&PRM_FLAG_SLIDER)
                {
                    float min = parameter.property(PP_MIN).toFloat();
                    float max = parameter.property(PP_MAX).toFloat();
                    QSlider *slider = (QSlider *)parameter.property(PP_WIDGET2).toLongLong();
                    float pos;
                    // value = min + pos/100(max-min)
                    //(value - min)100/(max-min) = pos
                    pos = (parameter.value().toFloat() - min)*SLIDER_SIZE/(max - min);
                    slider->setSliderPosition((int)pos);
                }
            }
        }
    }
    return 0;
}

void ConfigDialog::load()
{
    render(&m_interpreter->m_pixyParameters, NULL, m_pixyTabs);
    show();
}

void ConfigDialog::render(ParameterDB *data, QGridLayout *layout, QTabWidget *tabs)
{
    int i;
    QWidget *tab;

    DBG("rendering config...");
    Parameters &parameters = data->parameters();

    for (i=0; i<parameters.size(); i++)
    {
        Parameter &parameter = parameters[i];
        uint flags = parameter.property(PP_FLAGS).toInt();

        if (flags&PRM_FLAG_INTERNAL) // don't render!
            continue;

        PType type = parameter.type();
        QPushButton *button = NULL;
        QCheckBox *cbox = NULL;
        QSlider *slider = NULL;

        QLineEdit *line = new QLineEdit();
        QLabel *label = new QLabel(parameter.id());
        label->setToolTip(parameter.help());
        label->setAlignment(Qt::AlignRight);
        parameter.setProperty(PP_WIDGET, (qlonglong)line);
        line->setMinimumWidth(50);
        line->setMaximumWidth(75);
        line->setToolTip(parameter.help());

        if (type!=PT_INTS8) // make sure it's a scalar type
        {
            if (flags&PRM_FLAG_PATH)
            {
                button = new QPushButton("Change...");
                button->setProperty("Parameter", (qlonglong)&parameter);
                connect(button, SIGNAL(clicked()), this, SLOT(handleChangeClicked()));
                button->setToolTip("Select a new path");
                line->setMinimumWidth(200);
                line->setMaximumWidth(300);
                line->setText(parameter.value().toString());
            }
            else if (flags&PRM_FLAG_CHECKBOX)
            {
                cbox = new QCheckBox();
                cbox->setProperty("Parameter", (qlonglong)&parameter);
                cbox->setChecked(parameter.value().toBool());
                cbox->setToolTip(parameter.help());
                connect(cbox, SIGNAL(clicked()), this, SLOT(handleCheckBox()));
                parameter.setProperty(PP_WIDGET, (qlonglong)cbox);
                delete line;
                line = NULL;
            }
            else if (flags&PRM_FLAG_SLIDER)
            {
                float min = parameter.property(PP_MIN).toFloat();
                float max = parameter.property(PP_MAX).toFloat();
                slider = new QSlider(Qt::Horizontal);
                slider->setProperty("Parameter", (qlonglong)&parameter);
                slider->setMinimumWidth(SLIDER_SIZE);
                slider->setMaximumWidth(SLIDER_SIZE);
                slider->setRange(0, SLIDER_SIZE);
                slider->setSingleStep(1);
                slider->setToolTip(parameter.help());
                parameter.setProperty(PP_WIDGET2, (qlonglong)slider);
                float pos;
                // value = min + pos/100(max-min)
                //(value - min)100/(max-min) = pos
                pos = (parameter.value().toFloat() - min)*SLIDER_SIZE/(max - min);
                slider->setSliderPosition((int)pos);
                if (type==PT_FLT32)
                    line->setText(QString::number(parameter.value().toFloat(), 'f', 6));
                else
                    line->setText(QString::number(parameter.value().toInt()));
                connect(slider, SIGNAL(sliderMoved(int)), this, SLOT(handleSlider(int)));
            }
            else if (type==PT_STRING)
            {
                line->setMinimumWidth(200);
                line->setMaximumWidth(300);
                line->setText(parameter.value().toString());
            }
            else if (type==PT_FLT32)
            {
                float val = parameter.value().toFloat();
                line->setText(QString::number(val, 'f', 6));
            }
            else if (!(flags&PRM_FLAG_SIGNED))
            {
                uint val = parameter.value().toUInt();
                if (flags&PRM_FLAG_HEX_FORMAT)
                    line->setText("0x" + QString::number(val, 16));
                else
                    line->setText(QString::number(val));
            }
            else
            {
                int val = parameter.valueInt();
                line->setText(QString::number(val));
            }

            // deal with categories-- create category tab if needed
            if (tabs)
            {
                QString category = parameter.property(PP_CATEGORY).toString();
                if (category=="")
                    category = CD_GENERAL;
                tab = findCategory(category, tabs);
                if (tab==NULL)
                {
                    tab = new QWidget();
                    layout = new QGridLayout(tab);
                    tabs->addTab(tab, category);
                }
                else
                    layout = (QGridLayout *)tab->layout();
            }
            layout->addWidget(label, i, 0);
            if (cbox)
                layout->addWidget(cbox, i, 1);
            if (slider)
                layout->addWidget(slider, i, 2);
            if (line)
                layout->addWidget(line, i, 1);
            if (button)
                layout->addWidget(button, i, 2);
        }
    }

    if (tabs)
    {
        // set stretch on all tabs
        for (i=0; true; i++)
        {
            QWidget *tab = tabs->widget(i);
            if (tab)
            {
                ((QGridLayout *)tab->layout())->setRowStretch(100, 1);
                ((QGridLayout *)tab->layout())->setColumnStretch(100, 1);
            }
            else
                break;
        }
    }
    if (layout)
    {
        layout->setRowStretch(100, 1);
        layout->setColumnStretch(100, 1);
    }

    DBG("rendering config done");
}

void ConfigDialog::accept()
{
    if (updateDB()<0)
        return;
    m_interpreter->saveParams();
    QDialog::accept();
}

void ConfigDialog::reject()
{
    DBG("reject called");
    // clear all shadows
    m_interpreter->m_pixymonParameters->clearShadow();
    m_interpreter->m_pixyParameters.clearShadow();
    // at this point only shadow parameters that have been modified
    // have their dirty bits set
    // calling saveParames will save only shadowed parameters (we didn't call updatedb)
    // saving the parameters will cause the dirty bit to be set on pixy and all paramters will
    // be reloaded (remember that parameters on the pixy side don't have dirty bits, but instead
    // we have a global dirty bit.
    m_interpreter->saveParams();
    QDialog::reject();
}

void ConfigDialog::apply(QAbstractButton *button)
{
    if (button->text()=="Apply")
    {
        if (updateDB()<0)
            return;
        m_interpreter->saveParams();
    }
}

void ConfigDialog::handleChangeClicked()
{
    // get button that was clicked
    QPushButton *button = (QPushButton *)sender();
    // then grab parameter associated with the button
    Parameter *parameter = (Parameter *)button->property("Parameter").toLongLong();

    // bring up file dialog
    QFileDialog fd(this);
    fd.setWindowTitle("Please choose a folder");
    fd.setDirectory(parameter->value().toString());
    fd.setFileMode(QFileDialog::Directory);
    fd.setOption(QFileDialog::ShowDirsOnly);

    if (fd.exec())
    {
        QStringList dir = fd.selectedFiles();
        // update text box
        QLineEdit *line = (QLineEdit *)parameter->property(PP_WIDGET).toLongLong();
        line->setText(dir[0]);
    }
}


void ConfigDialog::handleCheckBox()
{
    QCheckBox *cbox = (QCheckBox *)sender();
    Parameter *parameter = (Parameter *)cbox->property("Parameter").toLongLong();

    m_interpreter->m_pixyParameters.mutex()->lock();
    parameter->set(cbox->isChecked(), true);
    parameter->setDirty(true);
    m_interpreter->m_pixyParameters.mutex()->unlock();
    m_interpreter->updateParam();
}

void ConfigDialog::handleSlider(int position)
{
    QSlider *slider = (QSlider *)sender();
    Parameter *parameter = (Parameter *)slider->property("Parameter").toLongLong();
    QLineEdit *line = (QLineEdit *)parameter->property(PP_WIDGET).toLongLong();
    float value;
    float min = parameter->property(PP_MIN).toFloat();
    float max = parameter->property(PP_MAX).toFloat();
    value = min + slider->sliderPosition()*(max - min)/SLIDER_SIZE;

    // use pixyParameters mutex as mutex between configdialog and rest of pixymon
    // namely worker thread in interpreter.  If we try to lock both mutexes
    // (pixymonParameters and pixyParameters) we can get into a double mutex deadlock
    // (as a rule, never lock more than 1 mutex at a time)
    m_interpreter->m_pixyParameters.mutex()->lock();
    if (parameter->type()==PT_FLT32)
    {
        parameter->set(value, true); // set as shadow
        line->setText(QString::number(value, 'f', 6));
    }
    else
    {
        parameter->set((int)value, true); // set as shadow
        line->setText(QString::number((int)value));
    }
    parameter->setDirty(true);
    m_interpreter->m_pixyParameters.mutex()->unlock();
    m_interpreter->updateParam();
}
