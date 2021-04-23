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
/*=========================================================================
 *
 * Eikonal Activation Simulation (EASI) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/


#ifndef QmitkCemrgRenderWindowEditor_h
#define QmitkCemrgRenderWindowEditor_h

#include <QmitkAbstractRenderEditor.h>

class QmitkRenderWindow;
class QmitkStdMultiWidget;

class QmitkCemrgRenderWindowEditor : public QmitkAbstractRenderEditor
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:

  berryObjectMacro(QmitkCemrgRenderWindowEditor)
  static const std::string EDITOR_ID;
  QmitkCemrgRenderWindowEditor();


  QmitkRenderWindow* GetActiveQmitkRenderWindow() const;
  QHash<QString,QmitkRenderWindow*> GetQmitkRenderWindows() const;
  QmitkRenderWindow* GetQmitkRenderWindow(const QString& id) const;
  mitk::Point3D GetSelectedPosition(const QString& id = QString()) const;
  void SetSelectedPosition(const mitk::Point3D &pos, const QString& id = QString());
  void EnableDecorations(bool enable, const QStringList& decorations = QStringList());
  bool IsDecorationEnabled(const QString& decoration) const;
  QStringList GetDecorations() const;

protected:
  virtual void CreateQtPartControl(QWidget *parent) override;
  virtual void SetFocus() override;

private:
  QmitkRenderWindow* m_RenderWindow0;
  QmitkRenderWindow* m_RenderWindow1;
  QmitkStdMultiWidget* m_standrdWindow;
};

#endif // QmitkCemrgRenderWindowEditor_h
