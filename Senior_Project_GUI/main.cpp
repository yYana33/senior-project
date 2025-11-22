#include <QtWidgets/QApplication>
#include "gui/mainwindow.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    MainWindow mainWindow;
    mainWindow.show();

    cout << "DNA Analyzer and Visualizer started successfully!" << endl;

    return app.exec();
}
