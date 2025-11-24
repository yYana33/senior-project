#include "MainWindow.h"
#include <iostream>
#include <qapplication.h>

using namespace std;

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setupUI();
    setupMenus();
    setupToolbar();

    connect(loadAction, &QAction::triggered, this, &MainWindow::loadFile);
    connect(loadSecondAction, &QAction::triggered, this, &MainWindow::loadSecondFile);
    connect(searchAction, &QAction::triggered, this, &MainWindow::searchSequence);
    connect(orfAction, &QAction::triggered, this, &MainWindow::findORFs);
    connect(alignAction, &QAction::triggered, this, &MainWindow::alignSequences);
    connect(aboutAction, &QAction::triggered, this, &MainWindow::showAbout);

    setWindowTitle("DNA Sequence Visualizer");
    setMinimumSize(1000, 700);

    sequenceDisplay->setPlaceholderText("No sequence loaded. Load a FASTA to begin.");
    resultsDisplay->setPlaceholderText("Analysis results will appear here.");
}

MainWindow::~MainWindow() {

}

void MainWindow::setupUI() {
    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    mainLayout = new QVBoxLayout(centralWidget);

    //the main content area
    contentLayout = new QHBoxLayout();

    sequenceDisplay = new QTextEdit(this);
    sequenceDisplay->setReadOnly(true);
    sequenceDisplay->setMaximumHeight(300);  //increased height so the info can fit
    sequenceDisplay->setFont(QFont("Courier", 9));

    resultsDisplay = new QTextEdit(this);
    resultsDisplay->setReadOnly(true);
    resultsDisplay->setFont(QFont("Courier", 9));
    resultsDisplay->setLineWrapMode(QTextEdit::NoWrap);

    contentLayout->addWidget(sequenceDisplay);
    contentLayout->addWidget(resultsDisplay);

    mainLayout->addLayout(contentLayout);

    statusBar()->showMessage("Ready to load DNA sequence");
}

void MainWindow::setupMenus() {
    QMenu* fileMenu = menuBar()->addMenu("&File");

    loadAction = new QAction("&Load FASTA...", this);
    loadAction->setShortcut(QKeySequence::Open);
    fileMenu->addAction(loadAction);

    loadSecondAction = new QAction("Load &Second FASTA...", this);
    loadSecondAction->setShortcut(Qt::CTRL | Qt::Key_2);
    fileMenu->addAction(loadSecondAction);

    QAction* testAction = new QAction("Load &Test Sequence", this);
    fileMenu->addAction(testAction);
    connect(testAction, &QAction::triggered, this, &MainWindow::loadTestSequence);

    fileMenu->addSeparator();

    exitAction = new QAction("E&xit", this);
    exitAction->setShortcut(QKeySequence::Quit);
    fileMenu->addAction(exitAction);
    connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

    fileMenu->addSeparator();

    exitAction = new QAction("E&xit", this);
    exitAction->setShortcut(QKeySequence::Quit);
    fileMenu->addAction(exitAction);
    connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

    QMenu* analysisMenu = menuBar()->addMenu("&Analysis");

    searchAction = new QAction("&Search Pattern...", this);
    searchAction->setShortcut(Qt::CTRL | Qt::Key_F);
    analysisMenu->addAction(searchAction);

    orfAction = new QAction("Find &ORFs", this);
    orfAction->setShortcut(Qt::CTRL | Qt::Key_G);
    analysisMenu->addAction(orfAction);

    alignAction = new QAction("&Align Sequences", this);
    alignAction->setShortcut(Qt::CTRL | Qt::Key_A);
    analysisMenu->addAction(alignAction);

    QMenu* helpMenu = menuBar()->addMenu("&Help");
    aboutAction = new QAction("&About", this);
    helpMenu->addAction(aboutAction);
}

void MainWindow::loadTestSequence() {
    updateDisplay();
    statusBar()->showMessage("Loaded test sequence for debugging");
}

void MainWindow::setupToolbar() {
    QToolBar* toolBar = addToolBar("Toolbar");

    toolBar->addAction(loadAction);
    toolBar->addAction(loadSecondAction);
    toolBar->addSeparator();
    toolBar->addAction(searchAction);
    toolBar->addAction(orfAction);
    toolBar->addAction(alignAction);
}

void MainWindow::loadFile() {
    QDir currentDir = QDir::current();
    qDebug() << "Current directory:" << currentDir.absolutePath();

    QString fileName = QFileDialog::getOpenFileName(
        this,
        "Open FASTA File",
        currentDir.absolutePath(),
        "FASTA Files (*.fasta *.fa *.fsa);;All Files (*.*)"
        );

    if (!fileName.isEmpty()) {
        statusBar()->showMessage("Loading: " + fileName);
        qDebug() << "User selected file:" << fileName;

        if (controller.loadSequence(fileName)) {
            statusBar()->showMessage("Successfully loaded: " + fileName);
            updateDisplay();
            qDebug() << "Displaying now";
        } else {
            statusBar()->showMessage("Failed to load: " + fileName);
            QMessageBox::warning(this, "Load error", "Could not load/validate the FASTA file.");
            qDebug() << "File load failed!";
        }
    } else {
        qDebug() << "User cancelled file selection";
    }
}

void MainWindow::loadSecondFile() {
    QDir currentDir = QDir::current();
    QString fileName = QFileDialog::getOpenFileName(
        this,
        "Open Second FASTA File for Alignment",
        currentDir.absolutePath(),
        "FASTA Files (*.fasta *.fa *.fsa);;All Files (*.*)"
        );

    if (!fileName.isEmpty()) {
        statusBar()->showMessage("Loading second sequence: " + fileName);

        if (controller.loadSecondSequence(fileName)) {
            statusBar()->showMessage("Successfully loaded second sequence: " + fileName);
            updateDisplay();
        } else {
            statusBar()->showMessage("Failed to load second sequence: " + fileName);
            QMessageBox::warning(this, "Load error", "Could not load/validate the second FASTA file.");
        }
    }
}

void MainWindow::updateDisplay() {
    qDebug() << "Updating display..";

    if (controller.hasSequence()) {
        QString displayText;

        //basic seq info
        displayText += controller.getSequenceInfo() + "\n\n";

        //stats
        displayText += controller.getStatisticsInfo() + "\n\n";

        //auto-search results
        if (!controller.hasSecondSequence()) {
            displayText += controller.getAutoSearchResults() + "\n\n";
        }

        sequenceDisplay->setPlainText(displayText);

        if (controller.hasSecondSequence()) {
            resultsDisplay->setPlainText(controller.getAlignmentInfo());
        } else {
            resultsDisplay->setPlainText(controller.getFeaturesInfo());
        }

    } else {
        sequenceDisplay->setPlainText("No sequence loaded. Use File → Load FASTA to begin.");
        resultsDisplay->setPlainText("Analysis results will appear here.");
    }

    sequenceDisplay->repaint();
    resultsDisplay->repaint();
}

void MainWindow::searchSequence() {
    if (!controller.hasSequence()) {
        QMessageBox::information(this, "No sequence", "Please load a sequence first.");
        return;
    }

    bool ok;
    QString pattern = QInputDialog::getText(this, "Search Pattern", "Enter DNA pattern to search for:", QLineEdit::Normal, "ATG", &ok);
    if (ok && !pattern.isEmpty()) {
        qDebug() << "MAINWINDOW starting search";
        controller.searchPattern(pattern);

        updateDisplay();

        QApplication::processEvents();

        statusBar()->showMessage("Search completed for: " + pattern);
        qDebug() << "MAINWINDOW search UI updated";
    }
}

void MainWindow::findORFs() {
    if (!controller.hasSequence()) {
        QMessageBox::information(this, "No sequence", "Please load a sequence first.");
        return;
    }

    qDebug() << "MAINWINDOW: Starting ORF finding";
    statusBar()->showMessage("Finding ORFs...");

    controller.findORFs();

    updateDisplay();

    QApplication::processEvents();

    statusBar()->showMessage("ORF finding completed");
    qDebug() << "MAINWINDOW: ORF UI updated";
}

void MainWindow::alignSequences() {
    if (!controller.hasSequence() || !controller.hasSecondSequence()) {
        QMessageBox::information(this, "Alignment",
        "Please load two sequences first:\n"
        "1. Use 'Load FASTA' for the first sequence\n"
        "2. Use 'Load Second FASTA' for the second sequence");
        return;
    }

    qDebug() << "MAINWINDOW starting sequence alignment";
    statusBar()->showMessage("Aligning sequences..");

    controller.alignSequences();

    updateDisplay();

    QApplication::processEvents();

    statusBar()->showMessage("Sequence alignment completed");
    qDebug() << "MAINWINDOW alignment UI updated";
}

void MainWindow::showAbout() {
    QMessageBox::about(this, "About DNA Visualizer",
                       "<h3>DNA Sequence Visualizer and Analyzer</h3>"
                       "<p>Senior Project by Yana Slavcheva</p>"
                       "<p>COS 4091A - Fall 2025</p>"
                       "<p>A tool for visualizing and analyzing DNA sequences with "
                       "Boyer-Moore search, ORF finding, and Smith-Waterman alignment.</p>");
}

