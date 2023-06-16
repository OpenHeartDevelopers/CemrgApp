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
CEMRG CMD APP TEMPLATE
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkCommandLineParser.h>
#include <mitkIOUtil.h>

// Qt
#include <QString>
#include <QDebug>
#include <QFileInfo>
#include <QProcess>

// C++ Standard
#include <algorithm>
#include <string>

// CemrgApp
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Tests");
    parser.setTitle("Template Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription(
        "This template command-line app for testing different functionalities.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    parser.addArgument(
        "reference", "i", mitkCommandLineParser::String,
        "Input object (reference)", "Reference for registration. Accepts any image format known to MITK.",
        us::Any(), false);
    parser.addArgument(
        "otherimage", "j", mitkCommandLineParser::String,
        "Input object (reference)", "Image to calculate transformation from. Accepts any image format known to MITK.",
        us::Any(), false);
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["reference"].Empty() ||
        parsedArgs["otherimage"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto input1 = us::any_cast<std::string>(parsedArgs["reference"]);
    auto input2 = us::any_cast<std::string>(parsedArgs["otherimage"]);

    // Default values for optional arguments
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose"))
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);

    try {
        // Code the functionality of the cmd app here.
        if (verbose) {
            MITK_INFO << "Verbose mode ON. Calculating registration.";
            MITK_INFO << "Input 1:" << input1;
            MITK_INFO << "Input 2:" << input2;
        }
        std::unique_ptr<QProcess> process(new QProcess());
        aPath = QCoreApplication::applicationDirPath() + "/MLib";
        process->setProcessChannelMode(QProcess::MergedChannels);
        process->setWorkingDirectory(aPath);

        QString INPUT1(input1.c_str());
        QString INPUT2(input2.c_str());

        //Setup registration
        QStringList arguments;
        QFileInfo fileinfo(INPUT1);

        QString dir = fileinfo.path();
        QString output = dir + "/rigid.dof";
        QString mirtk = aPath + "/register";

        MITK_INFO(verbose) << "OUTPUT VALUE: " << output.toStdString();

        arguments << INPUT1;
        arguments << INPUT2;
        arguments << "-dofout" << output;
        arguments << "-model" << "Rigid";
        arguments << "-verbose" << "3";

        MITK_INFO(verbose) << "Running command:" << mirtk.toStdString() << std::endl;
        process->start(mirtk, arguments);
        MITK_INFO(verbose) << "EXIT CODE: " << process->exitCode();
        MITK_INFO(verbose) << "Finished working. Output in file:\n\t" << output.toStdString() << std::endl;
        MITK_INFO(verbose) << "Closing Process" << std::endl;
        process->close();

    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
