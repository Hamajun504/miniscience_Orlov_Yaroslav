#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

double BACK_FRONT_Y_BOUND = 12.437;
double LEG_Z_BOUND = -61.765;
double FOOT_Z_BOUND = -70.36;
double LOWER_NECK_Z_BOUND = -15.308;
double UPPER_NECK_Z_BOUND = 19.752;
double MIDDLE = 12.266;
double LEG_AMP_MUL = 0.4;
double LEG_MEAN_LEN_MUL = 0.6;
double NECK_ANGLE_TANGENT = 0.4;
double CIRCLE_FREQ = 2 * M_PI;


enum class BodyPart {neck, head, foot_fl, foot_fr, foot_bl, foot_br, 
                leg_fl, leg_fr, leg_bl, leg_br, other};

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    double x0, y0, z0;

    BodyPart bodypart;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), x0(0), y0(0), z0(0), bodypart(BodyPart::other)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z) 
            : x(x), y(y), z(z), x0(x), y0(y), z0(z), bodypart(BodyPart::other)
    {
        define_bodypart();
        move(0);
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double time) {
        switch (bodypart)
        {
        case BodyPart::head:
        {
            x = x0 + (UPPER_NECK_Z_BOUND - LOWER_NECK_Z_BOUND) * NECK_ANGLE_TANGENT *
                cos(CIRCLE_FREQ * time);
            y = y0 + (UPPER_NECK_Z_BOUND - LOWER_NECK_Z_BOUND) * NECK_ANGLE_TANGENT *
                sin(CIRCLE_FREQ * time);
            break;
        }
        case BodyPart::neck:
        {
            x = x0 + (z0 - LOWER_NECK_Z_BOUND) * NECK_ANGLE_TANGENT *
                cos(CIRCLE_FREQ * time);
            y = y0 + (z0 - LOWER_NECK_Z_BOUND) * NECK_ANGLE_TANGENT *
                sin(CIRCLE_FREQ * time);
            break;
        }
        case BodyPart::foot_br:
        case BodyPart::foot_fl:
        {
            z = LEG_Z_BOUND + (FOOT_Z_BOUND - LEG_Z_BOUND) * (LEG_MEAN_LEN_MUL + LEG_AMP_MUL * cos(CIRCLE_FREQ * time)) + (z0 - FOOT_Z_BOUND);
            break;
        }
        case BodyPart::foot_bl:
        case BodyPart::foot_fr:
        {
            z = LEG_Z_BOUND + (FOOT_Z_BOUND - LEG_Z_BOUND) * (LEG_MEAN_LEN_MUL - LEG_AMP_MUL * cos(CIRCLE_FREQ * time)) + (z0 - FOOT_Z_BOUND);
            break;
        }
        case BodyPart::leg_br:
        case BodyPart::leg_fl:
        {
            z = LEG_Z_BOUND + (z0 - LEG_Z_BOUND) * (LEG_MEAN_LEN_MUL + LEG_AMP_MUL * cos(CIRCLE_FREQ * time));
            break;
        }
        case BodyPart::leg_bl:
        case BodyPart::leg_fr:
        {
            z = LEG_Z_BOUND + (z0 - LEG_Z_BOUND) * (LEG_MEAN_LEN_MUL - LEG_AMP_MUL * cos(CIRCLE_FREQ * time));
            break;
        }
        }
    }


private:
    void define_bodypart()
    {
        if (z > UPPER_NECK_Z_BOUND) bodypart = BodyPart::head;
        else if (z > LOWER_NECK_Z_BOUND) bodypart = BodyPart::neck;
        else if (z < FOOT_Z_BOUND)
        {
            if (y > BACK_FRONT_Y_BOUND & x > MIDDLE) bodypart = BodyPart::foot_bl;
            else if (y > BACK_FRONT_Y_BOUND & x < MIDDLE) bodypart = BodyPart::foot_br;
            else if (y < BACK_FRONT_Y_BOUND & x > MIDDLE) bodypart = BodyPart::foot_fl;
            else if (y < BACK_FRONT_Y_BOUND & x < MIDDLE) bodypart = BodyPart::foot_fr;
        }
        else if (z < LEG_Z_BOUND)
        {
            if (y > BACK_FRONT_Y_BOUND & x > MIDDLE) bodypart = BodyPart::leg_bl;
            else if (y > BACK_FRONT_Y_BOUND & x < MIDDLE) bodypart = BodyPart::leg_br;
            else if (y < BACK_FRONT_Y_BOUND & x > MIDDLE) bodypart = BodyPart::leg_fl;
            else if (y < BACK_FRONT_Y_BOUND & x < MIDDLE) bodypart = BodyPart::leg_fr;
        }
    }
};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;
    double time;

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints): time(0) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];

            nodes[i] = CalcNode(pointX, pointY, pointZ);

        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau) {
        time += tau;
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].move(time);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "movement_snaps/goose_movement_step_" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    // Шаг по времени
    double tau = 1. / 30;

    const unsigned int GMSH_TETR_CODE = 4;

    // Теперь придётся немного упороться:
    // (а) построением сетки средствами gmsh,
    // (б) извлечением данных этой сетки в свой код.
    gmsh::initialize();
    gmsh::open("goose.msh");

    //gmsh::fltk::run();

    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    // TODO: неплохо бы полноценно данные сетки проверять, да

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    for (size_t i = 0; i < 150; ++i)
    {
        mesh.doTimeStep(tau);
        mesh.snapshot(i);
        std::cout << "    \r";
        std::cout.flush();
        std::cout << round(i * 100. / 150) << "%\r";
        std::cout.flush();
    }

    return 0;
}
