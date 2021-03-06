#include "crpropa/ParticleState.h"

#include <HepPID/ParticleIDMethods.hh>

#include <stdexcept>
#include <algorithm>

namespace crpropa {

ParticleState::ParticleState() :
		id(0), charge(0), pmass(0), energy(0), position(0, 0, 0), direction(-1,
				0, 0) {
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getR();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy); // prevent negative energies
}

double ParticleState::getEnergy() const {
	return energy;
}

void ParticleState::setId(int newId) {
	id = newId;
	if (HepPID::isNucleus(id)) {
		pmass = nucleusMass(id);
		charge = HepPID::Z(id) * eplus; // HepPID::charge doesn't work for nuclei
		if (id < 0)
			charge *= -1; // HepPID::Z returns positive charge numbers for anti-nuclei
	} else {
		// pmass missing for non-nuclei
		charge = HepPID::charge(id) * eplus;
	}
}

int ParticleState::getId() const {
	return id;
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getCharge() const {
	return charge;
}

double ParticleState::getLorentzFactor() const {
	return energy / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	energy = lf * pmass * c_squared;
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3d ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

} // namespace crpropa
