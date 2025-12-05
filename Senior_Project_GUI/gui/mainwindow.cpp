#include "mainwindow.h"
#include <iostream>
#include <qapplication.h>
#include <QLabel>
#include <QRegularExpression>
#include <QActionGroup>

using namespace std;

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    sequenceVisualizer = new SequenceVisualizer(this);
    controller = new Controller(this);
    gcCircle = new GCContentCircle(this);

    setupUI();
    setupMenus();
    setupViewToolbar();

    connect(loadAction, &QAction::triggered, this, &MainWindow::loadFile);
    connect(loadSecondAction, &QAction::triggered, this, &MainWindow::loadSecondFile);
    connect(searchAction, &QAction::triggered, this, &MainWindow::searchSequence);
    connect(orfAction, &QAction::triggered, this, &MainWindow::findORFs);
    connect(alignAction, &QAction::triggered, this, &MainWindow::alignSequences);
    connect(aboutAction, &QAction::triggered, this, &MainWindow::showAbout);

    connect(clearSearchAction, &QAction::triggered, this, &MainWindow::clearSearchResults);
    connect(clearORFAction, &QAction::triggered, this, &MainWindow::clearORFResults);
    connect(clearAlignmentAction, &QAction::triggered, this, &MainWindow::clearAlignmentResults);
    connect(clearAllAction, &QAction::triggered, this, &MainWindow::clearAllResults);
    connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

    if (sequenceVisualizer) {
        connect(sequenceVisualizer, &SequenceVisualizer::viewModeChanged, this, &MainWindow::onViewModeChanged);
        connect(sequenceVisualizer, &SequenceVisualizer::zoomLevelChanged, this, &MainWindow::onZoomLevelChanged);
    }
    setWindowTitle("DNA Sequence Visualizer");
    setMinimumSize(1000, 700);

    resultsDisplay->setPlaceholderText("Load a sequence to view analysis.");

    updateDisplay();
}

MainWindow::~MainWindow() {

}
void MainWindow::setupUI() {
    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    //visualizer on top, text panel below
    mainLayout = new QVBoxLayout(centralWidget);
    mainLayout->setContentsMargins(8, 8, 8, 8);
    mainLayout->setSpacing(10);


    //visualization
    sequenceVisualizer->setMinimumHeight(180);
    sequenceVisualizer->setMaximumHeight(300);
    sequenceVisualizer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

    QHBoxLayout* topLayout = new QHBoxLayout();
    topLayout->addWidget(sequenceVisualizer, 3);
    topLayout->addWidget(gcCircle, 1);

    mainLayout->addLayout(topLayout);


    //results panel
    resultsDisplay = new QTextEdit(this);
    resultsDisplay->setReadOnly(true);
    resultsDisplay->setFont(QFont("Courier", 10));
    resultsDisplay->setLineWrapMode(QTextEdit::NoWrap);
    resultsDisplay->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    mainLayout->addWidget(resultsDisplay);

    statusBar()->showMessage("Ready to load DNA sequence");
}


void MainWindow::updateVisualization() {
    if (!sequenceVisualizer || !controller) return;

    sequenceVisualizer->setSequence(controller->getSequencePointer());

    sequenceVisualizer->setSearchMatches(
        controller->getSearchMatchPositions(),
        controller->getPatternLength()
        );

    std::vector<SequenceVisualizer::ORFRegion> visORFs;

    for (const auto& o : controller->getORFs()) {
        visORFs.push_back({ o.start, o.end, o.strand });
    }

    sequenceVisualizer->setORFs(visORFs);


    sequenceVisualizer->setAlignmentBlocks(controller->getAlignmentBlocks());

    bool hasSearch   = controller->hasSearchResults();
    bool hasORFs     = controller->hasORFResults();
    bool hasAlign    = controller->hasAlignmentResults();

    bool useAlignmentView = hasAlign && !hasSearch && !hasORFs;

    sequenceVisualizer->updateMode(
        controller->hasSearchResults(),
        controller->hasORFResults(),
        controller->hasAlignmentResults()
        );

    sequenceVisualizer->update();
    gcCircle->setGC(controller->getGCPercentage());

}

void MainWindow::setupMenus() {
    QMenu* fileMenu = menuBar()->addMenu("&File");

    //Loading operations
    loadAction = new QAction("&Load FASTA...", this);
    loadAction->setShortcut(QKeySequence::Open);
    loadAction->setStatusTip("Load a DNA sequence from FASTA file");
    fileMenu->addAction(loadAction);

    loadSecondAction = new QAction("Load &Second FASTA...", this);
    loadSecondAction->setShortcut(Qt::CTRL | Qt::Key_2);
    loadSecondAction->setStatusTip("Load a second sequence for alignment");
    fileMenu->addAction(loadSecondAction);

    fileMenu->addSeparator();

    //clear submenu
    QMenu* clearMenu = fileMenu->addMenu("&Clear Results");
    clearMenu->setStatusTip("Clear analysis results");

    clearSearchAction = new QAction("Clear &Search Results", this);
    clearSearchAction->setShortcut(Qt::CTRL | Qt::SHIFT | Qt::Key_S);
    clearSearchAction->setStatusTip("Remove all pattern search results");
    clearMenu->addAction(clearSearchAction);

    clearORFAction = new QAction("Clear &ORF Results", this);
    clearORFAction->setShortcut(Qt::CTRL | Qt::SHIFT | Qt::Key_O);
    clearORFAction->setStatusTip("Remove all ORF/gene predictions");
    clearMenu->addAction(clearORFAction);

    clearAlignmentAction = new QAction("Clear &Alignment Results", this);
    clearAlignmentAction->setShortcut(Qt::CTRL | Qt::SHIFT | Qt::Key_A);
    clearAlignmentAction->setStatusTip("Remove sequence alignment results");
    clearMenu->addAction(clearAlignmentAction);

    clearAllAction = new QAction("Clear &All Results", this);
    clearAllAction->setShortcut(Qt::CTRL | Qt::SHIFT | Qt::Key_C);
    clearAllAction->setStatusTip("Remove all analysis results");
    clearMenu->addAction(clearAllAction);

    fileMenu->addSeparator();

    //exit
    exitAction = new QAction("E&xit", this);
    exitAction->setShortcut(QKeySequence::Quit);
    exitAction->setStatusTip("Exit the application");
    fileMenu->addAction(exitAction);

    //analysis menu
    QMenu* analysisMenu = menuBar()->addMenu("&Analysis");
    analysisMenu->setStatusTip("Perform DNA sequence analysis");

    searchAction = new QAction("&Search Pattern...", this);
    searchAction->setShortcut(Qt::CTRL | Qt::Key_F);
    searchAction->setStatusTip("Search for a DNA pattern using Boyer-Moore algorithm");
    analysisMenu->addAction(searchAction);

    orfAction = new QAction("Find &ORFs", this);
    orfAction->setShortcut(Qt::CTRL | Qt::Key_G);
    orfAction->setStatusTip("Find Open Reading Frames (potential genes)");
    analysisMenu->addAction(orfAction);

    alignAction = new QAction("&Align Sequences", this);
    alignAction->setShortcut(Qt::CTRL | Qt::Key_A);
    alignAction->setStatusTip("Align two sequences using Smith-Waterman algorithm");
    analysisMenu->addAction(alignAction);

    QMenu* helpMenu = menuBar()->addMenu("&Help");

    aboutAction = new QAction("&About", this);
    aboutAction->setStatusTip("About DNA Sequence Visualizer");
    helpMenu->addAction(aboutAction);
}

void MainWindow::setupViewToolbar() {
    viewToolBar = addToolBar("View Controls");
    viewToolBar->setMovable(false);

    detailedViewAction = new QAction("Detailed", this);
    detailedViewAction->setCheckable(true);
    detailedViewAction->setToolTip("Show detailed view with colored bases");

    featureViewAction = new QAction("Feature", this);
    featureViewAction->setCheckable(true);
    featureViewAction->setToolTip("Show simplified feature view");

    autoViewAction = new QAction("Auto", this);
    autoViewAction->setCheckable(true);
    autoViewAction->setChecked(true);
    autoViewAction->setToolTip("Automatically switch between views based on zoom");

    //action group so only one can be selected at a time
    QActionGroup* viewGroup = new QActionGroup(this);
    viewGroup->addAction(detailedViewAction);
    viewGroup->addAction(featureViewAction);
    viewGroup->addAction(autoViewAction);

    zoomInAction = new QAction("Zoom In", this);
    zoomInAction->setShortcut(Qt::CTRL | Qt::Key_Plus);
    zoomInAction->setToolTip("Zoom in (Ctrl++)");

    zoomOutAction = new QAction("Zoom Out", this);
    zoomOutAction->setShortcut(Qt::CTRL | Qt::Key_Minus);
    zoomOutAction->setToolTip("Zoom out (Ctrl+-)");

    zoomResetAction = new QAction("Reset Zoom", this);
    zoomResetAction->setShortcut(Qt::CTRL | Qt::Key_0);
    zoomResetAction->setToolTip("Reset zoom to default (Ctrl+0)");

    //widgets to toolbar
    viewToolBar->addWidget(new QLabel("View Mode: "));
    viewToolBar->addAction(detailedViewAction);
    viewToolBar->addAction(featureViewAction);
    viewToolBar->addAction(autoViewAction);
    viewToolBar->addSeparator();
    viewToolBar->addAction(zoomInAction);
    viewToolBar->addAction(zoomOutAction);
    viewToolBar->addAction(zoomResetAction);

    //connect signals
    connect(detailedViewAction, &QAction::triggered, this, [&] {
        sequenceVisualizer->updateMode(false, false, false);
    });


    connect(featureViewAction, &QAction::triggered, this, [&] {
        sequenceVisualizer->updateMode(true, true, false);
    });


    connect(autoViewAction, &QAction::triggered, this, [&] {
        sequenceVisualizer->updateMode(
            controller->hasSearchResults(),
            controller->hasORFResults(),
            controller->hasAlignmentResults()
            );
    });

    connect(zoomInAction, &QAction::triggered, sequenceVisualizer, [this] {
        sequenceVisualizer->setZoomLevel(sequenceVisualizer->getZoomLevel() * 1.2);
    });


    connect(zoomOutAction, &QAction::triggered, this, [this]() {
        if (sequenceVisualizer) {
            sequenceVisualizer->setZoomLevel(sequenceVisualizer->getZoomLevel() / 1.5);
        }
    });

    connect(zoomResetAction, &QAction::triggered, this, [this]() {
        if (sequenceVisualizer) {
            sequenceVisualizer->setZoomLevel(1.0);
        }
    });
}

void MainWindow::onViewModeChanged(SequenceVisualizer::VisualizationMode mode) {

    detailedViewAction->setChecked(mode == SequenceVisualizer::DETAILED_VIEW);
    featureViewAction->setChecked(mode == SequenceVisualizer::FEATURE_VIEW);
    autoViewAction->setChecked(mode == SequenceVisualizer::ALIGNMENT_VIEW);

    QString modeText;

    switch (mode) {
    case SequenceVisualizer::DETAILED_VIEW:
        modeText = "Detailed View";
        break;

    case SequenceVisualizer::FEATURE_VIEW:
        modeText = "Feature View";
        break;

    case SequenceVisualizer::ALIGNMENT_VIEW:
        modeText = "Alignment View";
        break;
    }

    statusBar()->showMessage("View mode changed to: " + modeText);
}


void MainWindow::onZoomLevelChanged(double zoom) {
    statusBar()->showMessage(QString("Zoom level: %1x").arg(zoom, 0, 'f', 1));
}

void MainWindow::loadTestSequence() {
    updateDisplay();
    statusBar()->showMessage("Loaded test sequence for debugging");
}

void MainWindow::loadFile() {
    QDir currentDir = QDir::current();

    QString fileName = QFileDialog::getOpenFileName(
        this,
        "Open FASTA File",
        currentDir.absolutePath(),
        "FASTA Files (*.fasta *.fa *.fsa);;All Files (*.*)"
        );

    if (!fileName.isEmpty()) {
        statusBar()->showMessage("Loading: " + fileName);

        if (controller->loadSequence(fileName)) {
            statusBar()->showMessage("Successfully loaded: " + fileName);
            updateDisplay();
        } else {
            statusBar()->showMessage("Failed to load: " + fileName);
            QMessageBox::warning(this, "Load error", "Could not load/validate the FASTA file.");
        }
    } else {
        qDebug() << "User cancelled file selection";
    }
}

void MainWindow::loadSecondFile() {

    if (!controller->hasSequence()) {
        QMessageBox::warning(this, "Load Sequence First",
        "Please load the first DNA sequence before loading a second one for alignment.\n\n"
        "Use File → Load FASTA to load your first sequence.");
        return;
    }
    QDir currentDir = QDir::current();
    QString fileName = QFileDialog::getOpenFileName(
        this,
        "Open Second FASTA File for Alignment",
        currentDir.absolutePath(),
        "FASTA Files (*.fasta *.fa *.fsa);;All Files (*.*)"
        );

    if (!fileName.isEmpty()) {
        statusBar()->showMessage("Loading second sequence: " + fileName);

        if (controller->loadSecondSequence(fileName)) {
            statusBar()->showMessage("Successfully loaded second sequence: " + fileName);
            updateDisplay();
        } else {
            statusBar()->showMessage("Failed to load second sequence: " + fileName);
            QMessageBox::warning(this, "Load error", "Could not load/validate the second FASTA file.");
        }
    }
}

void MainWindow::updateDisplay() {
    if (!controller) return;

    if (controller->hasSequence()) {
        QString displayText;

        displayText += controller->getSequenceInfo() + "\n\n";

        //stats
        displayText += controller->getStatisticsInfo() + "\n\n";

        //auto-search results only when no second sequence yet
        if (!controller->hasSecondSequence()) {
            displayText += controller->getAutoSearchResults() + "\n\n";
        }

        //status
        displayText += "RESULTS\n";

        if (controller->hasSearchResults()) {
            displayText += QString("Search results: %1 matches\n").arg(controller->getSearchMatchCount());
        } else {
            displayText += "No search results\n";
        }

        if (controller->hasORFResults()) {
            displayText += QString("✓ ORF results: %1 genes found\n").arg(controller->getGeneCount());
        } else {
            displayText += "No ORF results\n";
        }

        if (controller->hasSecondSequence()) {
            if (controller->hasAlignmentResults()) {
                displayText += "Alignment results available\n";
            } else {
                displayText += "No alignment results (click Align)\n";
            }
        }

        if (controller->hasSecondSequence()) {
            displayText += "\n\nALIGNMENT DETAILS:\n";
            displayText += controller->getAlignmentInfo();
        } else {
            displayText += "\n\nFEATURE DETAILS:\n";
            displayText += controller->getFeaturesInfo();
        }

        resultsDisplay->setPlainText(displayText);

        sequenceVisualizer->updateMode(
            controller->hasSearchResults(),
            controller->hasORFResults(),
            controller->hasSecondSequence() && controller->hasAlignmentResults()
            );
        updateVisualization();
    } else {
        resultsDisplay->setPlainText(
            "No sequence loaded. Use File → Load FASTA to begin.\n\n"
            "Analysis results will appear here."
            );
        if (sequenceVisualizer) {
            sequenceVisualizer->clearVisualization();
        }
    }

    resultsDisplay->repaint();
}

void MainWindow::searchSequence() {
    if (!controller->hasSequence()) {
        QMessageBox::information(this, "No Sequence", "Please load a sequence first.");
        return;
    }

    //custom dialog with validation
    QDialog dialog(this);
    dialog.setWindowTitle("Search Pattern");

    QVBoxLayout* layout = new QVBoxLayout(&dialog);

    QLabel* label = new QLabel("Enter DNA pattern to search for:\n" "(Only A, T, C, G allowed, 2-20 characters)", &dialog);
    layout->addWidget(label);

    QLineEdit* patternEdit = new QLineEdit(&dialog);
    patternEdit->setPlaceholderText("e.g., ATG, TAA, GGCC");
    patternEdit->setMaxLength(20);
    layout->addWidget(patternEdit);

    QLabel* validationLabel = new QLabel(&dialog);
    validationLabel->setStyleSheet("color: red;");
    layout->addWidget(validationLabel);

    QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
    layout->addWidget(buttonBox);

    //connecting validation
    connect(patternEdit, &QLineEdit::textChanged, this, [patternEdit, validationLabel](const QString& text) {
        QString upperText = text.toUpper();

        //validating characters
        QRegularExpression dnaRegex("^[ATCG]*$");
        bool validChars = dnaRegex.match(upperText).hasMatch();

        //validating length
        bool validLength = text.length() >= 2 && text.length() <= 20;

        if (text.isEmpty()) {
            validationLabel->clear();
            patternEdit->setStyleSheet("");
        } else if (!validChars) {
            validationLabel->setText("Error! Only A, T, C, G characters allowed");
            patternEdit->setStyleSheet("border: 1px solid red;");
        } else if (!validLength) {
            validationLabel->setText("Error! Pattern must be 2-20 characters");
            patternEdit->setStyleSheet("border: 1px solid red;");
        } else {
            validationLabel->setText("Valid DNA pattern");
            validationLabel->setStyleSheet("color: green;");
            patternEdit->setStyleSheet("border: 1px solid green;");
        }
    });

    connect(buttonBox, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);

    //showing uppercase as the user types
    connect(patternEdit, &QLineEdit::textChanged, this, [patternEdit](const QString& text) {
        QSignalBlocker blocker(patternEdit);
        patternEdit->setText(text.toUpper());
    });

    if (dialog.exec() == QDialog::Accepted) {
        QString pattern = patternEdit->text().toUpper();


        if (pattern.length() < 2 || pattern.length() > 20) {
            QMessageBox::warning(this, "Invalid Pattern", "Pattern must be between 2 and 20 characters.");
            return;
        }

        //final validation
        QRegularExpression dnaRegex("^[ATCG]+$");
        if (!dnaRegex.match(pattern).hasMatch()) {
            QMessageBox::warning(this, "Invalid Pattern", "Pattern can only contain A, T, C, or G characters.");
            return;
        }

        controller->searchPattern(pattern);
        updateDisplay();
        updateVisualization();

        connect(controller, &Controller::errorOccurred, this, [this](const QString& error) {
            QMessageBox::warning(this, "Search Error", error);
        });

        controller->searchPattern(pattern);
        updateDisplay();

        //showing match count in status bar
        int matchCount = controller->getSearchMatchCount();

        statusBar()->showMessage(QString("Found %1 matches for: %2").arg(matchCount).arg(pattern));
    }
}

void MainWindow::findORFs() {
    if (!controller->hasSequence()) {
        QMessageBox::information(this, "No sequence", "Please load a sequence first.");
        return;
    }

    statusBar()->showMessage("Finding ORFs...");

    controller->findORFs();

    updateDisplay();

    QApplication::processEvents();

    statusBar()->showMessage("ORF finding completed");
}

void MainWindow::alignSequences() {
    if (!controller->hasSequence() || !controller->hasSecondSequence()) {
        QMessageBox::information(this, "Alignment",
        "Please load two sequences first:\n"
        "1. Use 'Load FASTA' for the first sequence\n"
        "2. Use 'Load Second FASTA' for the second sequence");
        return;
    }

    statusBar()->showMessage("Aligning sequences..");

    controller->alignSequences();

    updateDisplay();

    QApplication::processEvents();

    statusBar()->showMessage("Sequence alignment completed");
}

void MainWindow::clearSearchResults() {
    if (!controller->hasSequence()) {
        QMessageBox::information(this, "No Sequence", "No sequence loaded to clear results from.");
        return;
    }

    if (!controller->hasSearchResults()) {
        QMessageBox::information(this, "No Results", "No search results to clear.");
        return;
    }

    //simple confirmation for search results
    int result = QMessageBox::question(this, "Clear Search Results", "Are you sure you want to clear all search results?",
                                            QMessageBox::Yes | QMessageBox::No, QMessageBox::No);

    if (result == QMessageBox::Yes) {
        controller->clearSearchResults();
        updateDisplay();
        statusBar()->showMessage("Search results cleared");
    }
}

void MainWindow::clearORFResults() {
    if (!controller->hasSequence()) {
        QMessageBox::information(this, "No Sequence", "No sequence loaded to clear results from.");
        return;
    }

    if (!controller->hasORFResults()) {
        QMessageBox::information(this, "No Results", "No ORF results to clear.");
        return;
    }

    //simple confirmation for ORF results
    int result = QMessageBox::question(this, "Clear ORF Results", "Are you sure you want to clear all ORF (gene) results?",
                                       QMessageBox::Yes | QMessageBox::No, QMessageBox::No);

    if (result == QMessageBox::Yes) {
        controller->clearORFResults();
        updateDisplay();
        statusBar()->showMessage("ORF results cleared");
    }
}

void MainWindow::clearAlignmentResults() {
    if (!controller->hasSecondSequence()) {
        QMessageBox::information(this, "No Alignment", "No alignment results to clear.");
        return;
    }

    if (!controller->hasAlignmentResults()) {
        QMessageBox::information(this, "No Results", "No alignment results to clear.");
        return;
    }

    //simple confirmation for alignment results
    int result = QMessageBox::question(this, "Clear Alignment Results", "Are you sure you want to clear alignment results?",
                                       QMessageBox::Yes | QMessageBox::No, QMessageBox::No);

    if (result == QMessageBox::Yes) {
        controller->clearAlignmentResults();
        updateDisplay();
        statusBar()->showMessage("Alignment results cleared");
    }
}

void MainWindow::clearAllResults() {
    if (!controller->hasSequence() && !controller->hasAlignmentResults()) {
        QMessageBox::information(this, "No Results", "No results to clear.");
        return;
    }

    QString confirmation = "Are you sure you want to clear ";
    bool hasAnything = false;

    if (controller->hasSearchResults()) {
        confirmation += "search results";
        hasAnything = true;
    }

    if (controller->hasORFResults()) {
        if (hasAnything) confirmation += ", ";
        confirmation += "ORF results";
        hasAnything = true;
    }

    if (controller->hasAlignmentResults()) {
        if (hasAnything) confirmation += ", ";
        confirmation += "alignment results";
        hasAnything = true;
    }

    confirmation += "?";

    if (!hasAnything) {
        QMessageBox::information(this, "No Results", "No results to clear.");
        return;
    }

    //confirmation dialog
    int result = QMessageBox::question(this, "Clear All Results", confirmation + "\n\nThis action cannot be undone.",
                                       QMessageBox::Yes | QMessageBox::No, QMessageBox::No);

    if (result == QMessageBox::Yes) {
        controller->clearAllResults();
        updateDisplay();
        statusBar()->showMessage("All results cleared");
    }
}

void MainWindow::showAbout() {
    QMessageBox::about(this,
    "About DNA Visualizer",
    "<h3>DNA Sequence Visualizer and Analyzer</h3>"
    "<p>Senior Project by Yana Slavcheva</p>"
    "<p>COS 4091A - Fall 2025</p>"
    "<p>A tool for visualizing and analyzing DNA sequences with "
    "Boyer-Moore search, ORF finding, and Smith-Waterman alignment.</p>");
}

