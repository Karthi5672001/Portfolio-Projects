--Select * 
--From [Portfolio Project1]..coviddeaths$
--order by 3,4

--Select * 
--From [Portfolio Project1]..Covidvaccs$
--order by 3,4

---- Selecting the necessary Data

--Select Location,date, total_cases, new_cases, total_deaths, population
--From [Portfolio Project1]..coviddeaths$
--order by 1,2

----Looking at total_cases vs total_deaths in UK as percentage mortalility
--Select Location,date, total_cases, total_deaths, (total_deaths/total_cases)*100 as Mortalitypercent
--From [Portfolio Project1]..coviddeaths$
--Where Location like '%kingdom%'
--order by 1,2

---- total_cases vs population
--Select Location,date, total_cases, population, (total_cases/population)*100 as PopulationInfectionpercent
--From [Portfolio Project1]..coviddeaths$
--Where Location like '%kingdom%'
--order by 1,2

-- Looking at Countries with highest infection rate
Select Location, population, MAX(total_cases) as HighestInfectionCount,  MAX((total_cases/population))*100 as Casepercent
From [Portfolio Project1]..coviddeaths$
Group by location,population
order by Casepercent desc

-- Looking at Countries with highest death count
Select Location,  MAX(cast(total_deaths as int)) as HighestDeathCount
From [Portfolio Project1]..coviddeaths$
Where continent is not null
Group by location
order by HighestDeathCount desc

-- Looking at Continents with highest death count
Select location,  MAX(cast(total_deaths as int)) as HighestDeathCount
From [Portfolio Project1]..coviddeaths$
Where continent is null 
--Where Location is not like '%income%'
Group by location
order by HighestDeathCount desc

-- Global
--Global Mortalitypercent final
Select SUM(new_cases) as total_cases, SUM(cast(new_deaths as int)) as total_deaths, SUM(cast(new_deaths as int))/SUM(New_Cases)*100 as Mortalitypercent
From [Portfolio Project1]..coviddeaths$
Where continent is not null
--Group By date
order by 1,2

--Global Mortalitypercent over time
Select date, SUM(new_cases) as total_cases, SUM(cast(new_deaths as int)) as total_deaths, SUM(cast(new_deaths as int))/SUM(New_Cases)*100 as Mortalitypercent
From [Portfolio Project1]..coviddeaths$
Where continent is not null
Group By date
order by 1,2

--Total population vs Vaccinations
Select dea.continent, dea.location, dea.date, population, vac.new_vaccinations
, SUM(CONVERT(int, vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.Date ) as Sumofpeoplevaccinated
From [Portfolio Project1]..Covidvaccs$ vac
join [Portfolio Project1]..coviddeaths$ dea
on dea.location = vac.location
and dea.date = vac.date
Where dea.continent is not null
order by 2,3

--Use CTE

with PopsvsVacs (Continent, Location, Date, Population, New_Vaccinations, Sumofpeoplevaccinated) as
(

Select dea.continent, dea.location, dea.date, population, vac.new_vaccinations
, SUM(CONVERT(int, vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.Date ) as Sumofpeoplevaccinated
From [Portfolio Project1]..Covidvaccs$ vac
join [Portfolio Project1]..coviddeaths$ dea
on dea.location = vac.location
and dea.date = vac.date
Where dea.continent is not null
--order by 2,3

)
Select *, (Sumofpeoplevaccinated/Population)*100 as Sumofpeoplevaccinatedpercent
From PopsvsVacs

--Temp table

DROP Table if exists #Sumofpeoplevaccinatedpercent
Create Table #Sumofpeoplevaccinatedpercent
(
Continent nvarchar(255),
Location nvarchar(255),
Date datetime,
Population numeric,
new_vaccinations numeric,
Sumofpeoplevaccinated numeric
)
Insert into #Sumofpeoplevaccinatedpercent
Select dea.continent, dea.location, dea.date, population, vac.new_vaccinations
, SUM(CONVERT(int, vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.Date ) as Sumofpeoplevaccinated
From [Portfolio Project1]..Covidvaccs$ vac
join [Portfolio Project1]..coviddeaths$ dea
on dea.location = vac.location
and dea.date = vac.date
Where dea.continent is not null
--order by 2,3

Select *, (Sumofpeoplevaccinated/Population)*100 as Sumofpeoplevaccinatedpercent
From #Sumofpeoplevaccinatedpercent

--Creating View to store data for visulisations

Create View Sumofpeoplevaccinatedpercent as
Select dea.continent, dea.location, dea.date, population, vac.new_vaccinations
, SUM(CONVERT(int, vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.Date ) as Sumofpeoplevaccinated
From [Portfolio Project1]..Covidvaccs$ vac
join [Portfolio Project1]..coviddeaths$ dea
on dea.location = vac.location
and dea.date = vac.date
Where dea.continent is not null
--order by 2,3