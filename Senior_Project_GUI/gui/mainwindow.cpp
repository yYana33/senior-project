#include "MainWindow.h"
#include <iostream>

using namespace std;

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setupUI();
    setupMenus();
    setupToolbar();

    connect(loadAction, &QAction::triggered, this, &MainWindow::loadFile);
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

    contentLayout = new QHBoxLayout();//main area

    //will upgrade
    sequenceDisplay = new QTextEdit(this);
    sequenceDisplay->setReadOnly(true);
    sequenceDisplay->setMaximumHeight(200);

    resultsDisplay = new QTextEdit(this);
    resultsDisplay->setReadOnly(true);

    contentLayout->addWidget(sequenceDisplay);
    contentLayout->addWidget(resultsDisplay);

    mainLayout->addLayout(contentLayout);

    statusBar()->showMessage("Ready to load a DNA sequence");
}

void MainWindow::setupMenus() {
    QMenu* fileMenu = menuBar()->addMenu("&File");

    loadAction = new QAction("&Load FASTA...", this);
    loadAction->setShortcut(QKeySequence::Open);
    fileMenu->addAction(loadAction);

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

    alignAction = new QAction("&Align Sequences...", this);
    alignAction->setShortcut(Qt::CTRL | Qt::Key_A);
    analysisMenu->addAction(alignAction);

    QMenu* helpMenu = menuBar()->addMenu("&Help");
    aboutAction = new QAction("&About", this);
    helpMenu->addAction(aboutAction);
}

void MainWindow::setupToolbar() {
    QToolBar* toolBar = addToolBar("Toolbar");

    toolBar->addAction(loadAction);
    toolBar->addSeparator();
    toolBar->addAction(searchAction);
    toolBar->addAction(orfAction);
    toolBar->addAction(alignAction);
}

void MainWindow::loadFile() {
    QString fileName = QFileDialog::getOpenFileName(
        this,
        "Open FASTA File",
        "",
        "FASTA Files (*.fasta *.fa *.fsa);;All Files (*.*)"
        );

    if (!fileName.isEmpty()) {
        statusBar()->showMessage("Loading: " + fileName);

        //TEMPORARY
        sequenceDisplay->setPlainText("Loaded: " + fileName + "\n\n[FASTA content would be displayed here]");
        resultsDisplay->setPlainText("File loaded successfully!\n\nUse Analysis menu to perform operations.");

        cout << "Loading file: " << fileName.toStdString() << endl;
    }
}

void MainWindow::searchSequence() {
    bool ok;
    QString pattern = QInputDialog::getText(this, "Search Pattern", "Enter DNA pattern to search for:", QLineEdit::Normal, "ATG", &ok);
    if (ok && !pattern.isEmpty()) {
        resultsDisplay->append("Searching for pattern: " + pattern);
        //Will update so it connects to the actual search algorithm
    }
}

void MainWindow::findORFs() {
    resultsDisplay->append("Finding ORFs...");
    //Will update
}

void MainWindow::alignSequences() {
    resultsDisplay->append("Aligning sequences...");
    //Will update
}

void MainWindow::showAbout() {
    QMessageBox::about(this, "About DNA Visualizer",
                       "<h3>DNA Sequence Visualizer and Analyzer</h3>"
                       "<p>Senior Project by Yana Slavcheva</p>"
                       "<p>COS 4091A - Fall 2025</p>"
                       "<p>A tool for visualizing and analyzing DNA sequences with "
                       "Boyer-Moore search, ORF finding, and Smith-Waterman alignment.</p>");
}
