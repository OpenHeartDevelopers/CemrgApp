/*=========================================================================
  
Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef QmitkCemrgRenderWindowEditor_H_
#define QmitkCemrgRenderWindowEditor_H_

#include <QmitkAbstractRenderEditor.h>

class QmitkRenderWindow;
class QmitkStdMultiWidget;


class QmitkCemrgRenderWindowEditor : public QmitkAbstractRenderEditor
{
    Q_OBJECT

public:

    berryObjectMacro(QmitkCemrgRenderWindowEditor)

    static const std::string EDITOR_ID;

    QmitkCemrgRenderWindowEditor();


    // -------------------  mitk::IRenderWindowPart  ----------------------

    /**
   * \see mitk::IRenderWindowPart::GetActiveQmitkRenderWindow()
   */
    QmitkRenderWindow* GetActiveQmitkRenderWindow() const;

    /**
   * \see mitk::IRenderWindowPart::GetQmitkRenderWindows()
   */
    QHash<QString,QmitkRenderWindow*> GetQmitkRenderWindows() const;

    /**
   * \see mitk::IRenderWindowPart::GetQmitkRenderWindow(QString)
   */
    QmitkRenderWindow* GetQmitkRenderWindow(const QString& id) const;

    /**
   * \see mitk::IRenderWindowPart::GetSelectionPosition()
   */
    mitk::Point3D GetSelectedPosition(const QString& id = QString()) const;

    /**
   * \see mitk::IRenderWindowPart::SetSelectedPosition()
   */
    void SetSelectedPosition(const mitk::Point3D &pos, const QString& id = QString());

    /**
   * \see mitk::IRenderWindowPart::EnableDecorations()
   */
    void EnableDecorations(bool enable, const QStringList& decorations = QStringList());

    /**
   * \see mitk::IRenderWindowPart::IsDecorationEnabled()
   */
    bool IsDecorationEnabled(const QString& decoration) const;

    /**
   * \see mitk::IRenderWindowPart::GetDecorations()
   */
    QStringList GetDecorations() const;

protected:

    void SetFocus();

    void CreateQtPartControl(QWidget* parent);

private:

    QmitkRenderWindow* m_RenderWindow0;
    QmitkRenderWindow* m_RenderWindow1;
    QmitkStdMultiWidget* m_standrdWindow;

};

#endif /*QmitkCemrgRenderWindowEditor_H_*/
